<?php
/*Check that this file is being accessed by the template*/
if (!isset($in_template))
  {
    header( 'Location: /index.php/404');
    return;
  }

$pagetitle="Documentation";
ob_start();
   ?>
<h1>Tutorials</h1>

<h1>Source Code</h1>

<!-- Page End -->
<?php $content = ob_get_clean(); ?>
