<?php
/*Check that this file is being accessed by the template*/
if (!isset($in_template))
  {
    header( 'Location: /index.php/404');
    return;
  }

$pagetitle="Tutorial 1: Compiling and Installing DynamO";
ob_start();
   ?>
<!-- Page Begin -->
<p>This tutorial covers the requirements, compilation and installation
of the DynamO simulation package. It is recommended that you build
your own version of DynamO to keep up with the rapid code development
and to ensure compatibility with your system.</p>
<h2>Step 0: Build Requirements</h2>
<p>Currently DynamO will only run on <strong>Gnu/Linux</strong> based
systems (e.g., Ubuntu/Gentoo/RedHat). You will need to be familiar
with how to install programs on whichever distribution of Linux you
are using before you will be able to setup DynamO. </p>
<p>DynamO, like many Linux programs, is driven through a Command-Line
Interface (CLI). To be able to use DynamO, you will need to be
familiar with the terminal of your Linux distribution. Take a look
at <a href="https://help.ubuntu.com/community/UsingTheTerminal">this
link</a> to learn more about the terminal and how it works if you are
at all unsure what this means.</p>
<p>Before you can build DynamO, you will need several programs and
libraries installed. There are also several optional libraries which,
if they're installed, will activate extra features. </p>
<h3>Required Libraries</h3>
<ul>
  <li>
    <a href="http://www.bzip.org/">libbz2</a>
    and <a href="http://www.zlib.net/">libz</a> - The output of DynamO
    is compressed for efficiency using these libraries.
    <br/>(<strong>Ubuntu Packages</strong>: libbz2-dev
    zlib1g-dev).<br />You will need the static version of these two
    libraries (this can be a problem on Suse linux).
  </li>
</ul>
<h3>Visualisation Requirements (Optional)</h3>
<ul>
    <li><a href="http://www.gtkmm.org/">Gtkmm</a>
    <br/>(<strong style="font-size: 16px;">Ubuntu Package</strong>:
    libgtkmm-2.4-dev).
    </li>
    <li>
      <a href="http://freeglut.sourceforge.net/">Freeglut</a>
      <br/>(<strong>Ubuntu Package</strong>:freeglut3-dev)
    </li>
    <li>
      <a href="http://glew.sourceforge.net/">GLEW</a> (ver 1.6 and
      above) <br/> (<strong>Ubuntu Package</strong>: libglew1.6-dev)
    </li>
    <li>
      <a href="http://www.khronos.org/opencl/">OpenCL</a> - An OpenCL
      implementation is provided with the latest AMD and NVidia
      graphics card drivers. You will need a relatively modern
      graphics card to use the visualiser too. <br/> (<strong>Ubuntu
      Packages</strong>: either fglrx-dev (AMD) or nvidia-dev
      (NVidia)).
    </li>
    <li>
      <a href="http://ffmpeg.org/">libavcodec</a> Allows you to record visualisations directly to a movie file.
      <br/> (<strong>Ubuntu Package</strong>: libavcodec-dev)
    </li>
</ul>
<h2>Step 1: Download the Source Code</h2>
<p> Use the menu link to the left to download a copy of DynamO. Once
you have the source code, change into the directory ready to start the
build.</p>
<div class="code">cd DynamO</div>
<h2>Step 2: Compilation and Installation</h2>
<p>
  DynamO uses the modern, powerful, but quite complicated boost-build
  system. Using the boost build system takes some getting used to;
  however, to make it easy to build DynamO there is a fake Makefile
  included in the sources.
</p>
<p>
  Building DynamO is then as straightforward as running the make
  command:
</p>
<div class="code">make</div>
<p>
  This step can take a while, it will download a copy of boost, and
  build DynamO.
</p>
<p>
  If there are any errors, they are often due to missing build
  dependencies. Dynamo automatically checks if it can find the
  dependencies it needs to build. The list of tests should look like
  this:
</p>
<div class="code">
  <p>Performing configuration checks</p>
  <p><strong>..The dependencies below are required to build DynamO..</strong><br />
    - DynamO: Static bzip2 library : yes<br /> - DynamO: Static zlib library : yes</p>
  <p><strong>..The dependencies below are for coil/visualisation
    only..</strong><br /> - Coil:
    Gtkmm&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
    : yes<br /> - Coil: OpenCL
    lib&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; : yes<br /> -
    Coil:
    GLEW&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
    : yes<br /> - Coil:
    GLUT&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
    : yes
  <p>- DynamO-Coil Integration&nbsp; : yes</p>
  <p><strong>..The dependencies below add extra functionality to the
    visualizer and are optional..</strong><br /> - Coil: libCwiid
    Wii-mote support (Optional) : yes<br /> - Magnet: libavcodec
    (video encoding support, may also be packaged with FFMPEG) :
    yes</p></p>
</div>
<p>
  If you are missing static versions of <strong>bzip2</strong>
  or <strong>zlib</strong>, then DynamO won't build at all. If you are
  missing any of Coil's dependencies [DynamO-Coil Integration
  : <strong>no</strong>] DynamO will still build, but without the
  visualizer support. 
</p>
<p>
  If you still have errors, take a look at the
  <a href="/index.php/documentation">documentation</a> to find ways of
  contacting the developers.
</p>
<h2>Step 3: Installing the Executables</h2>
<p>
  Once DynamO compiled successfully the two main exectuables,
  called <em><strong>dynamod</strong></em>
  and <em><strong>dynarun</strong></em>, should be in the <em>bin</em>
  directory. You can run these programs from there, or copy them to
  some convenient place. For example, you can look at the help page of
  dynamod by running the following command from the <em>DynamO</em> directory.:
</p>
<div class="code">./bin/dynamod --help</div>
<p>
   You can also install the exectuables into your /usr/bin directory
  (although its not recommended) using the following command.
</p>
<div class="code">sudo make install</div>
<p>Congratulations! You now have a working installation of DynamO.</p>
<h2>Step 4: Updating the Code</h2>
<p>
  You can easily update to the latest version of DynamO, just enter
  the DynamO directory again and run the following command
</p>
<div class="code">git pull</div>
<p>
  Then just run <em>"make"</em> again:</p>
  <div class="code">make</div> 
<p>(and "<em>sudo make install"</em> if you previously installed it) to build
new executables with the latest changes.
</p>
<h2>Appendix A: Building Executables for Debugging</h2>
<p>If you're having some trouble with DynamO, you can build a debug
version of the simulator once a normal version has been built. This
version is slower than the standard version, but contains many extra
sanity checks and verbose error reports. To create the debugging
version just run the following command in the DynamO directory</p>
<div class="code">src/boost/bjam -j2 install debug</div>
<p>
  This will install some executables built with debugging symbols and
  extra sanity checks in the <em>bin/</em> directory. These
  executables have the suffix "<em>_d</em>" (dynamod_d and dynarun_d)
  to indicate they're the debugging version.
</p>
<!-- Page End -->
<?php $content = ob_get_clean(); ?>
