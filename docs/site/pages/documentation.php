<?php
/*Check that this file is being accessed by the template*/
if (!isset($in_template))
  {
    header( 'Location: /index.php/404');
    return;
  }

$pagetitle="Documentation";
ob_start();
   ?>
<h1>Tutorials</h1> 
<p>
  The documentation for DynamO is still being written and limited to a
  set of tutorials on the basic topics. Please click on any of the
  links below to take a look
</p>
<div style="text-align:center;"><?php button("Tutorial 1: Compiling DynamO from Source","/index.php/tutorial1");?></div>
<p>
  If what you need to simulate is not covered in the tutorials, but is
  listed in the <a href="/index.php/features">features</a>, please
  feel free to email the developers for some advice (see below).
</p>

<h1>Source Code</h1>
<p>
  If you're looking to extend DynamO or to understand how it works,
  you'll need to take a look at the source code. The DynamO API is
  partially documented using Doxygen and a up to date version is
  available at the link below.
</p>
<div style="text-align:center;"><?php button("DynamO API Documentation","/doxygen");?></div>
</table> 

<h1>Support</h1> 

<p>If you have a problem and cannot find the answer in the
documentation, you can email your queries to the following
address:</p>
<p style="text-align:center;">support@dynamomd.org</p>
<!-- Page End -->
<?php $content = ob_get_clean(); ?>
