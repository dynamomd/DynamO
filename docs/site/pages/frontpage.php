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
<p>It can be used to simulate a wide range of particle models, from hard
spheres, to square-wells and stepped potentials, thin needles and
more. From molecular dynamics to granular dynamics, DynamO attempts to
implement the latest in event-driven algorithms. Take a look at the
<a href="/index.php/features">features</a> to see if DynamO does what you need.</p>

<!-- Page End -->
<h1>Latest News</h1>
<?php 
 //Look in the news directory and create a date sorted list of the news items
 $file_array=array();
 foreach (glob('pages/news/*.html') as $filename)
 {
    $file_array[filectime($filename)]=basename($filename); // or just $filename
 }
 ksort($file_array);

 if (!empty($file_array))
 { include('pages/news/'.current($file_array)); }

 $content = ob_get_clean(); 
?>
