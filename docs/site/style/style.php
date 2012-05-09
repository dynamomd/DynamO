<?php header("Content-type: text/css"); 
$footerheight = "31px";
?>

/* Main html element styles */
html { height: 100%; }
/*#e4f4ff*/
body { margin: 0; padding: 0; height: 100%; background: #000000; }

/* Floating footer code */
#wrapper { min-height:100%; overflow-y:visible; }
#wrapperfooterpad { height:<?php echo $footerheight; ?>; width:100%; margin-bottom:10px; clear:both; }

#header { position:relative; padding-left:30px; padding-right:30px; height:115px;}
#sitelogo { background-repeat:no-repeat; background-image:url(../images/sitelogo.png); background-position: 50px 0; }

#footer {
    height:<?php echo $footerheight; ?>; 
    margin-top:-<?php echo $footerheight; ?>;
    position:relative;
}

.w3footerlogo { border:0; width:88px; height:31px; position:absolute; }

.sprite { background-image: url(../images/csssprites.png); }

/*     Styling of the round edged boxes       */
.borderleft { left:15px; position:absolute; top:0px; width:15px; height:100%; }
.borderright { right:15px; position:absolute; top:0px; width:15px; height:100%; }
.bordercentre { background-color:#ffffff; position:absolute; left:30px; right:30px; bottom:0; top:0;}
.topleftcornerborder { width:15px; height:15px; position:absolute; top:0; background-position: 0px 0px; }
.bottomleftcornerborder { width:15px; height:15px; position:absolute; bottom:0px; background-position: 0px -15px; }
.toprightcornerborder { width:15px; height:15px; position:absolute; top:0; background-position: -15px 0px; }
.bottomrightcornerborder { width:15px; height:15px; position:absolute; bottom:0px; background-position: -15px -15px; }
.verticalborder { background: #ffffff; position:absolute; top:15px; bottom:15px;left:0px; right:0px; }