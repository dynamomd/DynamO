<?php header("Content-type: text/css"); 
$footerheight = "35px";
?>

/* Main html element styles */
html { height: 100%; }
body { margin: 0; padding: 0; height: 100%; background: #e7f4fe; }

/* Floating footer code */
#wrapper { min-height:100%; overflow-y:visible; }
#wrapperfooterpad { height:<?php echo $footerheight; ?>; width:100%; clear:both; }

#leftmenu { width:20%; min-height:100%; float:left; display:block; }
#content { width:80%; display:block; float:left; }

#footer { 
    text-align:right;
    height:<?php echo $footerheight; ?>; 
    margin-top:-<?php echo $footerheight; ?>;
/*    background: url("../images/footerback.png") repeat-x scroll 0 0 transparent;*/
}

.footerlogo { display: block; float: right; border:0; width:88px; height:31px }