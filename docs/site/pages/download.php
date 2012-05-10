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
<h1>Download</h1>

<a class="button" href="https://github.com/toastedcrumpets/DynamO">
  <div class="left"></div><div class="center">Download</div><div class="right"></div>
</a>
<!--- Page End --->
<?php $content = ob_get_clean(); ?>
