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
Test content<br/>More test content<br/>
Test content<br/>More test content<br/>
Test content<br/>More test content<br/>
Test content<br/>More test content<br/>
Test content<br/>More test content<br/>
<!--- Page End --->

<?php $content = ob_get_clean(); ?>
