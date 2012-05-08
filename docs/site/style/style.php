<?php header("Content-type: text/css"); 
$footerheight = "35px";
?>

/* Main html element styles */
html { height: 100%; }
/*#e4f4ff*/
body { margin: 0; padding: 0; height: 100%; background: #000000; }

/* Floating footer code */
#wrapper { min-height:100%; overflow-y:visible; }
#wrapperfooterpad { height:<?php echo $footerheight; ?>; width:100%; clear:both; }

#leftmenu {  }
#content { }

#footer {
    height:<?php echo $footerheight; ?>; 
    margin-top:-<?php echo $footerheight; ?>;
    position:relative;
}

.footerlogo { border:0; width:88px; height:31px }
.sprite { background-image: url(../images/csssprites.png); }