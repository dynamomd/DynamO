<?php
/*Check that this file is being accessed by the template*/
if (!isset($in_template))
  {
    header( 'Location: /index.php/404');
    return;
  }

$pagetitle="Features";
ob_start();
   ?>
<!--- Page Begin --->
<div class="jce_caption" style="width: 365px; float: right; display: inline-block;"><iframe src="http://www.youtube.com/embed/VGSY6DNrrxk" width="360" height="220" marginheight="0" marginwidth="0" scrolling="auto" frameborder="0" align="bottom" name="plg_iframe">No Frames</iframe>
<div style="clear: both;">A polydisperse, sheared system undergoing compression.</div>
</div>
<h2>Simple Systems</h2>
<p>DynamO can simulate:</p>
<ul>
<li><strong>Smooth</strong> or <strong>rough hard spheres</strong>, <strong>square wells</strong>, <strong>soft cores</strong> and all the other simple EDMD potentials.</li>
<li><strong>Millions of particles </strong>for <strong>billions of events</strong>: Using 500 bytes per particle and running at around 75k events a second, even an old laptop can run huge simulations with accurate results.</li>

<li><strong>Non-Equilibrium Molecular Dynamics</strong> (NEMD): DynamO has Lees-Edwards boundary conditions for shearing systems, thermalised walls and Andersen or rescaling thermostats.</li>
<li><strong>Compression dynamics</strong>: Need high density systems? DynamO can compress any of its systems using isotropic compaction until the pressure diverges! </li>
<li><strong>Poly-dispersity</strong>: All interactions are generalised to fully poly-disperse particles.</li>
</ul>

<div class="jce_caption" style="width: 365px; float: right; display: inline-block;"><iframe src="http://www.youtube.com/embed/9oaobaxhGX8" width="360" height="220" marginheight="0" marginwidth="0" scrolling="auto" frameborder="0" align="bottom" name="plg_iframe">No Frames</iframe>
<div style="clear: both;">Hopper emptying through a chute modelled with sphere meshes, and accelerated using sleeping particles (dynamic particles turn white when sleeping).</div>
</div>

<h2>Granular Systems/Complex Boundaries</h2>
<p>DynamO contains all the latest EDMD algorithms for dissipative systems:</p>
<ul>
<li><strong>Event driven dynamics with gravity</strong>: Almost all interactions work in the presence of gravity, allowing event driven simulations of macroscopic systems.</li>
<li><strong>Sleeping particles</strong>: In systems where particles come to a rest, the sleeping particles algorithm can remove the cost of simulating the particles completely.</li>
<li><strong>Triangle or sphere meshes</strong>: Arbitrary complex boundaries can be implemented using triangle meshes imported from CAD drawings.</li>
<li><strong>Stepped potentials</strong>: Arbitrary stepped potentials may be simulated to approximate any soft potential.</li>
</ul>

<h2>Polymeric Systems/Accelerated Dynamics</h2>
<p>Dynamo has extremely advanced algorithms for accelerated dynamics:</p>
<ul>
<li><strong>Parallel tempering/Replica exchange</strong>: Run multiple simulations in parallel and use this Monte Carlo technique to enhance the sampling of phase space.</li>
<li><strong>Histogram reweighting</strong>: Combine the results from replica exchange simulations and extrapolate to temperatures which were not simulated.</li>
<li><strong>Multicanonical simulations</strong>: Used in conjunction with replica exchange techniques, multicanonical simulations greatly enhance the sampling of phase space, helping find the true native state of the polymer.</li>
</ul>

<!--- Page End --->
<?php $content = ob_get_clean(); ?>
