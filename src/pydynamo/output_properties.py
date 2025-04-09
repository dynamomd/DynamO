import os
import pickle
from subprocess import check_call

import numpy
import scipy

from pydynamo.config_files import ConfigFile
from pydynamo.file_types import XMLFile, validate_xmlfile
from pydynamo.weighted_types import WeightedArray, WeightedFloat


# A XMLFile/ElementTree but specialised for DynamO output files
class OutputFile(XMLFile):
    def __init__(self, filename):
        super().__init__(filename)
        
    def N(self):
        return int(self.tree.find('.//ParticleCount').attrib['val'])

    def events(self):
        return int(self.tree.find('.//Duration').attrib['Events'])

    def t(self):
        return float(self.tree.find('.//Duration').attrib['Time'])

    def numdensity(self):
        return float(self.tree.find('.//Density').attrib['val'])

    def __getitem__(self, key):
        return OutputFile.output_props[key](self)

    def __str__(self):
        return "OutputFile("+self._filename+")"
    
    output_props = {}

def validate_outputfile(filename):
    return validate_xmlfile(filename)
    
class OutputProperty:
    def __init__(self, dependent_statevars : list, dependent_outputs : list, dependent_outputplugins: list):
        self._dep_statevars = dependent_statevars
        self._dep_outputs = dependent_outputs
        self._dep_outputplugins = dependent_outputplugins

    def init(self):
        return None

    def result(self, state, outputfile, configfilename, counter, manager, output_dir):
        return None
    
class SingleAttrib(OutputProperty):
    def __init__(self, tag, attrib, dependent_statevars, dependent_outputs, dependent_outputplugins, time_weighted=True, div_by_N=False, div_by_t=False, missing_val = 0, skip_missing=False):
        OutputProperty.__init__(self, dependent_statevars, dependent_outputs, dependent_outputplugins)
        self._tag = tag
        self._attrib = attrib
        self._time_weighted = time_weighted
        self._div_by_N = div_by_N
        self._div_by_t = div_by_t
        self._missing_val = missing_val
        self._skip_missing=skip_missing

    def init(self):
        return WeightedFloat()

    def value(self, outputfile):
        tag = outputfile.tree.find('.//'+self._tag)
        if tag is None:
            if self._missing_val is None:
                raise RuntimeError('Failed to find the tag "'+self._tag+'" in the outputfile')
            else:
                return self._missing_val
            
        if self._attrib not in tag.attrib:
            if self._missing_val is None:
                raise RuntimeError('Failed to find attribute "'+self._attrib+'" in the tag "'+self._tag+'" in the outputfile')
            else:
                return self._missing_val
                
        val = float(tag.attrib[self._attrib])
        if self._div_by_N:
            N = outputfile.N()
            if N > 0:
                val /= outputfile.N()
        if self._div_by_t:
            t = outputfile.t()
            if t > 0:
                val /= t
        return val

    def weight(self, outputfile):
        tag = outputfile.tree.find('.//'+self._tag)
        #If we have a missing_val defined, then we can use and weight it, else don't give this any weight
        if ((self._missing_val is None) or self._skip_missing) and (tag is None or self._attrib not in tag.attrib):
            return 0
        
        if self._time_weighted:
            return float(outputfile.tree.find('.//Duration').attrib['Time'])
        else:
            return float(outputfile.tree.find('.//Duration').attrib['Events'])

    def result(self, state, outputfile, configfilename, counter, manager, output_dir):
        return WeightedFloat(self.value(outputfile), self.weight(outputfile))

def parseToArray(text):
    data = []
    for row in text.split('\n'):
        row_data = list(map(float, row.split()))
        if len(row_data) > 0:
            data.append(row_data)
    return numpy.array(data)

class CollisionMatrixOutputProperty(OutputProperty):
    def __init__(self):
        OutputProperty.__init__(self, dependent_statevars=[], dependent_outputs=[], dependent_outputplugins=['-LCollisionMatrix'])

    def result(self, state, outputfile, configfilename, counter, manager, output_dir):
        return None

class VACFOutputProperty(OutputProperty):
    def __init__(self):
        OutputProperty.__init__(self, dependent_statevars=[], dependent_outputs=[], dependent_outputplugins=['-LVACF'])

    def result(self, state, outputfile, configfilename, counter, manager, output_dir):
        restart_idx = output_dir.split('_')[-1]
        filename_root = manager.workdir+'_VACF/' + manager.statename(state, var_separator='/') + "/run_" + restart_idx + '_' + str(counter)
        os.makedirs(filename_root, exist_ok=True)
        for tag in outputfile.tree.findall('.//VACF/Particles/Species'):
            pickle.dump(parseToArray(tag.text), open(filename_root + '/species_'+tag.attrib['Name']+'.pkl', 'wb'))
        
        for tag in outputfile.tree.findall('.//VACF/Topology/Structure'):
            pickle.dump(parseToArray(tag.text), open(filename_root + '/topology_'+tag.attrib['Name']+'.pkl', 'wb'))
        return None

class RadialDistributionOutputProperty(OutputProperty):
    def __init__(self):
        OutputProperty.__init__(self, dependent_statevars=[], dependent_outputs=[], dependent_outputplugins=['-LRadialDistribution'])

    def init(self):
        return WeightedArray()

    def result(self, state, outputfile, configfilename, counter, manager, output_dir):
        #Presume that each tag is in order, and has a common bin width
        samples = float(outputfile.tree.find('.//RadialDistributionMoments').attrib["SampleCount"])
        bin_width = 0.01        
        #Grab all the moments, drop the R values for now
        moment_tags = outputfile.tree.findall('.//RadialDistributionMoments/Species/Moment')
        moments = numpy.array([[float(line.split()[1]) for line in tag.text.strip().split("\n")] for tag in moment_tags])

        # In the simulation we have collected the moments of, N(r), the number of pairs below a radius r. 
        # These moments are collected about an origin/offset N₀(r), i.e. ⟨(N(r)-N₀(r))^n⟩
        # The origin is a numerical trick to reduce the size of the moments to reduce precision issues, its an
        # approximation of the mean, by taking the initial value of N(r) at t=0.

        #moments is a 6xNbins array with N0, then <(N-N0)>, <(N-N0)^2>, and <(N-N0)^3> etc.
        N0 = numpy.copy(moments[0,:])
        #We set the N0 column (copied above) into the <(N-N0)^0> column. This makes the formula exact 
        moments[0,:] = numpy.ones_like(moments[0,:])

        # We eventually want to convert to the central moments ⟨(N(r)-⟨N(r)⟩)^n⟩, To
        # do this we need to determine the first moment ⟨N(r)⟩. There is a general
        # expression for transforming moments:
        #
        # ⟨(x-b)^n⟩ = Σ_{i=0}^n comb(n, i) ⟨(x-a)^i⟩ (a-b)^{n-i}
        #
        # comb is combination function in scipy.special (binomial coefficient).
        # But the answer for this first step is simple, x=N(r), b=0, a=N₀(r), n=1 gives the straightforward identity
        #
        # ⟨N(r)⟩ = ⟨N(r)-N₀(r)⟩ + N₀(r)

        # For all other steps n>2, b=⟨N(r)⟩, and a = N₀(r).
        # ⟨(N(r)-⟨N(r)⟩)^n⟩ = Σ_{i=0}^n comb(n, i) ⟨(N(r)-N₀(r))^i⟩ (N₀(r)-⟨N(r)⟩)^{n-i}

        central_moments = numpy.zeros((moments.shape[0]-1, moments.shape[1]))
        central_moments[0,:] = moments[1,:] + N0
        ##Now make the other central moments
        for n in range(2, moments.shape[0]):
            central_moments[n-1] = sum(scipy.special.comb(n, i) * moments[i] * (N0-central_moments[0,:])**(n-i) for i in range(n+1))

        return WeightedArray(central_moments, samples)

class RadialDistEndOutputProperty(OutputProperty):
    def __init__(self):
        OutputProperty.__init__(self, dependent_statevars=[], dependent_outputs=[], dependent_outputplugins=[])


    def result(self, state, outputfile, configfilename, counter, manager, output_dir):
        #This ending is not used, we only want one output for speed.
        #restart_idx = output_dir.split('_')[-1]
        #+ "/run_" + restart_idx + '_' + str(counter)
        filename_root = manager.workdir+'_RadialDist/' + manager.statename(state, var_separator='/')

        if os.path.isdir(filename_root):
            return None

        os.makedirs(filename_root, exist_ok=True)
        
        logfile = open(os.path.join(filename_root, 'run.log'), 'a')
        check_call(["dynarun", configfilename, '-o', os.path.join(filename_root, 'RadDist.config.xml.bz2'), '-c', '0', "--out-data-file", os.path.join(filename_root, 'RadDist.out.xml.bz2'), '-LRadialDistribution'], stdout=logfile, stderr=logfile)
        
        of = OutputFile(os.path.join(filename_root, 'RadDist.out.xml.bz2'))
        for tag in of.tree.findall('.//RadialDistribution/Species'):
            A = tag.attrib['Name1']
            B = tag.attrib['Name2']
            output_pkl = filename_root + '/species_'+A+'_'+B+'.pkl'
            pickle.dump(parseToArray(tag.text), open(output_pkl, 'wb'))
        return None

class OrderParameterProperty(OutputProperty):
    '''
    See here https://freud.readthedocs.io/en/stable/modules/order.html#freud.order.Steinhardt
    
    This property is expensive to run at data collection time, as it processes every configuration file to determine the order parameter.
    The advantage is that it can be run on any simulation, no need for extra output plugins.
    '''
    def __init__(self, L):
        OutputProperty.__init__(self, dependent_statevars=[], dependent_outputs=[], dependent_outputplugins=[])
        self.L = L
    
    def result(self, state, outputfile, configfilename, counter, manager, output_dir):
        # We use freud to calculate the Steinhardt order parameter
        # This function is called on every "run" of the simulation, with the outputfile and config file
        # So we process the particle positions in the config file to calculate the order parameter
        import freud
        configfile = ConfigFile(configfilename)
        box, points = configfile.to_freud()

        #Steinhardt for FCC
        ql = freud.order.Steinhardt(self.L)
        ql.compute((box, points), neighbors={"num_neighbors": self.L})
        ql_value = ql.particle_order
        
        return WeightedFloat(numpy.mean(ql_value), 1)

    def init(self):
        return WeightedFloat()
    
    def result(self, state, outputfile, configfilename, counter, manager, output_dir):
        import freud
        configfile = ConfigFile(configfilename)
        box, points = configfile.to_freud()

        #Steinhardt for FCC
        ql = freud.order.Steinhardt(self.L)
        ql.compute((box, points), neighbors={"num_neighbors": self.L})
        ql_value = ql.particle_order
        
        return WeightedFloat(numpy.mean(ql_value), 1)

class ChungLuConfigurationModel(OutputProperty):
    '''
    https://en.wikipedia.org/wiki/Configuration_model#Chung-Lu_Configuration_Model
    

    '''
    def __init__(self, L):
        OutputProperty.__init__(self, dependent_statevars=[], dependent_outputs=[], dependent_outputplugins=[])
        self.L = L
    
    def result(self, state, outputfile, configfilename, counter, manager, output_dir):
        
        configfile = ConfigFile(configfilename)
        if len(configfile.tree.findall('.//Interaction/CaptureMap')) != 1:
            raise RuntimeError("ChungLuConfigurationModel requires a single CaptureMap in the config file")

        import networkx as nx
        G = nx.Graph()
        for particle in configfile.tree.findall('.//Pt'):
            G.add_node(int(particle.attrib['ID']))

        for pair in configfile.tree.findall('.//Interaction/CaptureMap/Pair'):
            G.add_edge(int(pair.attrib['ID1']), int(pair.attrib['ID2']))
        
        return WeightedFloat(nx.community.modularity(G), 1)

    def init(self):
        return WeightedFloat()
    
    def result(self, state, outputfile, configfilename, counter, manager, output_dir):
        import freud
        configfile = ConfigFile(configfilename)
        box, points = configfile.to_freud()

        #Steinhardt for FCC
        ql = freud.order.Steinhardt(self.L)
        ql.compute((box, points), neighbors={"num_neighbors": self.L})
        ql_value = ql.particle_order
        
        return WeightedFloat(numpy.mean(ql_value), 1)


OutputFile.output_props["N"] = SingleAttrib('ParticleCount', 'val', [], [], [], missing_val=None)#We use missing_val=None to cause an error if the tag is missing
OutputFile.output_props["p"] = SingleAttrib('Pressure', 'Avg', [], [], [], missing_val=None)
OutputFile.output_props["cv"] = SingleAttrib('ResidualHeatCapacity', 'Value', [], [], [], div_by_N=True, missing_val=None)
OutputFile.output_props["u"] = SingleAttrib('UConfigurational', 'Mean', [], [], [], div_by_N=True, missing_val=None)
OutputFile.output_props["T"] = SingleAttrib('Temperature', 'Mean', [], [], [], missing_val=None)
OutputFile.output_props["density"] = SingleAttrib('Density', 'val', [], [], [], missing_val=None)
OutputFile.output_props["MSD"] = SingleAttrib('MSD/Species', 'diffusionCoeff', [], [], ['-LMSD'], missing_val=None, skip_missing=True)
OutputFile.output_props["NeventsSO"] = SingleAttrib('EventCounters/Entry[@Name="SOCells"]', # Outputfile tag name
                                                    'Count', # Outputfile tag attribute name
                                                    ["Rso"], # Required state variable
                                                    [], # Required output variables
                                                    [], # Required output plugins
                                                    div_by_N=True, # Divide the count by N
                                                    div_by_t=True, # Also divide by t
                                                    missing_val=0) # If counter is missing, return 0
OutputFile.output_props["VACF"] = VACFOutputProperty()
OutputFile.output_props["RadialDistEnd"] = RadialDistEndOutputProperty()
OutputFile.output_props["RadialDistribution"] = RadialDistributionOutputProperty()
OutputFile.output_props["FCCOrder"] = OrderParameterProperty(6)
OutputFile.output_props["CollisionMatrix"] = CollisionMatrixOutputProperty()
