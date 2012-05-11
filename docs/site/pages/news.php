<?php
/*Check that this file is being accessed by the template*/
if (!isset($in_template))
  {
    header( 'Location: /index.php/404');
    return;
  }

$pagetitle="News";

//Look in the news directory and create a date sorted list of the news items
$file_array=glob('pages/news/*.html');
rsort($file_array);
ob_start();
foreach ($file_array as $filename) { echo "<hr />"; include($filename); }
$content = ob_get_clean();
?>
