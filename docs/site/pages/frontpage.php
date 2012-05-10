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

<!--- Page Begin --->
<!--- Page End --->

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
