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

#header { position:relative; height:115px; }

#footer {
    height:<?php echo $footerheight; ?>; 
    margin-top:-<?php echo $footerheight; ?>;
    position:relative;
}

.footerlogo { border:0; width:88px; height:31px }

.sprite { background-image: url(../images/csssprites.png); }

/*     Styling of the round edged boxes       */
.borderleft { left:15px; position:absolute; top:0px; width:15px; height:100%; }
.borderright { right:15px; position:absolute; top:0px; width:15px; height:100%; }
.bordercentre { height:100%; background-color:#ffffff; position:absolute; left:30px; right:30px;}
.topleftcornerborder { width:15px; height:15px; position:absolute; top:0; background-position: 0px 0px; }
.bottomleftcornerborder { width:15px; height:15px; position:absolute; bottom:0px; background-position: 0px -15px; }
.toprightcornerborder { width:15px; height:15px; position:absolute; top:0; background-position: -15px 0px; }
.bottomrightcornerborder { width:15px; height:15px; position:absolute; bottom:0px; background-position: -15px -15px; }
.verticalborder { background: #ffffff; position:absolute; top:15px; bottom:15px;left:0px; right:0px; }