<?php
/*Check that this file is being accessed by the template*/
if (!isset($in_template))
  {
    header( 'Location: /index.php/404');
    return;
  }

$pagetitle="Welcome";
$content="Test content<br/>More test content<br/>
	  Test content<br/>More test content<br/>
	  Test content<br/>More test content<br/>
	  Test content<br/>More test content<br/>
	  Test content<br/>More test content<br/>
"
?>