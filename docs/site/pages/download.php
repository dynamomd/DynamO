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
<p>DynamO is an open-source and free code, made available under
version 3 of the GPL licence. There are several ways you can obtain a
copy of the source code and these are listed below.</p>

<h1>Git Access to the Source</h1> 

<p>The recommended method to downloade an up-to-date copy of the
DynamO sources is using git and the <a
href="https://github.com/toastedcrumpets/DynamO"> public GitHub
repository</a>.</p>

<p> Just use the following git command in your terminal:</p> 

<div class="code"> git clone
https://github.com/toastedcrumpets/DynamO.git
</div>

<p>This will create a folder called <em>DynamO</em> in the working
directory.</p>
<p>
  You can now pick which branch or version of DynamO you'd like. There
  are several available, but the default <strong>master</strong>
  branch is probably the most popular. There are three main branches
  available:
</p>
<ul>
    <li>
      <strong>master:</strong> This is the stable development branch,
      code should compile and run without errors but new features may
      be a little buggy.
    </li>
    <li>
      <strong>dynamo-1-2:</strong> This is the most recent stable
      release, the code is as stable as possible but may be up three
      months behind in features.
    </li>
    <li>
      <strong>experimental:</strong> This is where new patches and
      features are born. Probably a bit unstable and buggy, but very
      fresh! Not really meant for the public, but its been made
      available to let others track or join in the development
      process.
    </li>
</ul>
<p>If you decide you want something other than the default branch
(master), just check it out using the following command:</p>
<div class="code">git checkout dynamo-1-2</div>

<h1>Alternative Source Code Downloads</h1>
<p>
  It is recommended that you use the git to download the source code
  whenever possible. However, incase there is any problem accessing
  git on your computer, you may download a copy of the sources using
  the links below.
</p>

<div style="text-align:center; padding: 5px;"><?php button("Download Latest Development Branch","https://github.com/toastedcrumpets/DynamO/zipball/master");?></div>
<div style="text-align:center; padding: 5px;"><?php button("Download Stable Version 1.3","https://github.com/toastedcrumpets/DynamO/zipball/dynamo-1-3");?></div>
<div style="text-align:center; padding: 5px;"><?php button("Download Stable Version 1.2","https://github.com/toastedcrumpets/DynamO/zipball/dynamo-1-2");?></div>

<!-- Page End -->
<?php $content = ob_get_clean(); ?>
