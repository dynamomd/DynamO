<?php
/*Check that this file is being accessed by the template*/
if (!isset($in_template))
  {
    header( 'Location: /index.php/404');
    return;
  }

$pagetitle="Download";
ob_start();
   ?>
<!-- Page Begin -->
<h1>Ubuntu/Gentoo packages</h1>
There are no live packages built yet for the latest version of the code.

<h1>Git Access to the Source</h1>
<p>You can download an up-to-date copy of the DynamO sources from
the<a href="https://github.com/toastedcrumpets/DynamO"> public GitHub
repository</a> using git. Just use the following git command:</p>
<div class="codemain">git clone https://github.com/toastedcrumpets/DynamO.git</div>
<p>This will create a folder called <em>DynamO</em> in the working
directory. If you have trouble using Git, you can directly download a
zip file of the latest sources using the links below.
<h1>Alternative Source Code Downloads</h1> 

<p>
  It is recommended that you use the git to download the source code
when possible. However, incase there is any problem accessing git on
your computer, you may download a copy of the sources using the links
below.</p>
<table border="0" padding="1">
  <tr>
    <td>
      <a class="button" href="https://github.com/toastedcrumpets/DynamO/zipball/master">
	<div class="left"></div><div class="center">Download</div><div class="right"></div>
      </a>
    </td>
    <td>Latest Development Branch</td>
  </tr>
  <tr>
    <td>
      <a class="button" href="https://github.com/toastedcrumpets/DynamO/zipball/dynamo-1-3">
	<div class="left"></div><div class="center">Download</div><div class="right"></div>
      </a>
    </td>
    <td>Version 1.3 Stable Code</td>
  </tr>
  <tr>
    <td>
      <a class="button" href="https://github.com/toastedcrumpets/DynamO/zipball/dynamo-1-2">
	<div class="left"></div><div class="center">Download</div><div class="right"></div>
      </a>
    </td>
    <td>Version 1.2 Stable Code</td>
  </tr>
</table>

<!-- Page End -->
<?php $content = ob_get_clean(); ?>
