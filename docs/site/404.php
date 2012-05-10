<?php
/*Check that this file is being accessed by the template*/
if (!isset($in_template))
  {
    header( 'Location: /index.php/404');
    return;
  }

$pagetitle="404 Error";
ob_start();
   ?>
<!--- Page Begin --->
Could not find the page you were looking for!
<!--- Page End --->
<?php $content = ob_get_clean(); ?>
