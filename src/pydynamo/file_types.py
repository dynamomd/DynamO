import numpy

try:
    import lxml.etree as ET

    #print("Running with lxml.etree")
except ImportError:
    try:
        # normal cElementTree install
        import xml.etree.cElementTree as ET

        #print("Running with cElementTree")
    except ImportError:
        try:
            # normal ElementTree install
            import xml.etree.ElementTree as ET

            import elementtree.ElementTree as ET

            #print("Running with ElementTree")
        except ImportError:
            print("Failed to import ElementTree from any known place")

class XMLFile:
    """A wrapper around  to allow loading and saving to
    bzip2 compressed files. """
    
    def __init__(self, filename, compressed=None):
        """Loads the xml file, decompressing it with bz2 first if the
        filename ends with .bz2"""
        import io
        if isinstance(filename, str):
            self._filename = filename        
            if  filename.endswith('.xml.bz2'):
                import bz2
                f = bz2.BZ2File(filename)
                self.tree = ET.parse(f)
                f.close()
            elif filename.endswith('.xml'):
                self.tree = ET.parse(filename)
            else:
                raise RuntimeError('Unknown file extension for configuration file load "'+filename+'"')
        elif isinstance(filename, io.BytesIO):
            data = filename.read()
            try:
                self._filename = filename.name
            except:
                self._filename = "bytestream.xml"
            try:
                import bz2
                data = bz2.decompress(data)
                print("Successfully decompressed!")
            except:
                print("Could not decompress, trying direct reading")
            self.tree = ET.fromstring(data)
        else:
            raise RuntimeError("Could not determine file type")
        
    def save(self, filename):
        if filename.endswith('.xml.bz2'):
            import bz2
            f = bz2.BZ2File(filename, mode='w')
            f.write(ET.tostring(self.tree.getroot()))
            f.close()
        elif filename.endswith('.xml'):
            open(filename, 'w').write(ET.tostring(self.tree.getroot()))
        else:
            raise RuntimeError('Unknown file extension for configuration file save "'+filename+'"')

    def __str__(self):
        return "XMLFile("+self._filename+")"

def validate_xmlfile(filename):
    import bz2
    try:
        ET.parse(bz2.BZ2File(filename))
        return True
    except Exception as e:
        print("#!!!#", filename, e)
        return False
    