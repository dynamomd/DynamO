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
<!--- Page Begin --->
<h1 style="text-align:center;">Download</h1>
<!--- Page End --->
<?php $content = ob_get_clean(); ?>
