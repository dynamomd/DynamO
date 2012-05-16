<?php
/*Check that this file is being accessed by the template*/
if (!isset($in_template))
  {
    header( 'Location: /index.php/404');
    return;
  }

$pagetitle="Welcome";
ob_start();
   ?>

<!-- Page Begin -->

<p>DynamO is a free, open-source, event-driven particle simulator.</p>
<p>It can be used to simulate a wide range of particle models, from
hard spheres, to square-wells and stepped potentials, thin needles and
more. From molecular dynamics to granular dynamics, DynamO implements
the latest in event-driven algorithms. Take a look at the <a
href="/index.php/features">features</a> to see if DynamO does what you
need.</p>

<!-- Page End -->
<h1>Latest News</h1>
<?php 
 //Look in the news directory and create a date sorted list of the news items
 $file_array=glob('pages/news/*.html');
 rsort($file_array);

 if (!empty($file_array))
 { include(current($file_array)); }
 $content = ob_get_clean(); 
?>
