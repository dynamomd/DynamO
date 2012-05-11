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
<table border="0" cellpadding="5px" style="margin-left:auto; margin-right:auto; text-align:center;">
  <tr>
    <td>
      <a class="button" href="/index.php/tutorial1">
	<div class="left"></div><div class="center">Tutorial 1: Compiling DynamO from Source</div><div class="right"></div>
      </a>
    </td>
  </tr>
</table> 
<p>
  If what you need to simulate is not covered in the tutorials, but is
  listed in the <a href="/index.php/features">features</a>, please
  feel free to ask the developers for some advice.
</p>

<h1>Source Code</h1>
<p>
  If you're looking to extend DynamO or to understand how it works,
  you'll need to take a look at the source code. The DynamO API is
  partially documented using Doxygen and a up to date version is
  available at the link below.
</p>
<table border="0" cellpadding="5px" style="margin-left:auto; margin-right:auto; text-align:center;">
  <tr>
    <td>
      <a class="button" href="/doxygen">
	<div class="left"></div><div class="center">DynamO API Documentation</div><div class="right"></div>
      </a>
    </td>
  </tr>
</table> 
<p>
  If what you need to simulate is not covered in the tutorials, but is
  listed in the <a href="/index.php/features">features</a>, please
  feel free to ask the developers for some advice.
</p>


<h1>Support</h1>

<!-- Page End -->
<?php $content = ob_get_clean(); ?>
