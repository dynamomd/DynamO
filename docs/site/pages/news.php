<?php
/*Check that this file is being accessed by the template*/
if (!isset($in_template))
  {
    header( 'Location: /index.php/404');
    return;
  }

$pagetitle="News";

//Look in the news directory and create a date sorted list of the news items
$file_array=array();
foreach (glob('pages/news/*.html') as $filename)
{
   $file_array[filectime($filename)]=basename($filename); // or just $filename
}
krsort($file_array);

//Now output the news items in order
ob_start();
foreach ($file_array as $filename)
{
   include('pages/news/'.$filename);
}
$content = ob_get_clean();
?>
