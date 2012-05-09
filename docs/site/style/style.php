<?php header("Content-type: text/css"); 
$pagebg = "#b0d0ff";
?>

/* Main html element styles */
html { height: 100%; }
/*#e4f4ff*/
body { margin: 0; padding: 0; height: 100%; background: <?php echo $pagebg; ?>; }

/* Floating footer code */
#wrapper { min-height:100%; overflow-y:visible; position:relative; }

#wrapperfooterpad { height:55px; width:100%; margin-bottom:10px; clear:both; }

#header { position:relative; padding-left:30px; padding-right:30px; height:145px;}
#sitelogo { background-repeat:no-repeat; background-image:url(../images/sitelogo.png); background-position: 15px 15px; }

#headercontentspacing { height: 15px; clear: both; }

#footer {
    height:40px; 
    margin-top:-55px;
    position:relative;
}

.w3footerlogo { border:0; width:88px; height:31px; position:absolute; top:5px; }
.sprite { background-image: url(../images/csssprites.png); } 

#leftmenu { width:200px;float:left; margin-left:15px; position:relative; }
#leftmenu a { 
    line-height:29px; 
    color:#000000; 
    text-decoration: none; 
    display:block; 
    padding-left:16px; 
    background-color:#ffffff; 
}
#leftmenu a:hover { background-color:#eeeeee; }


#contentwrapper { margin-left:230px; margin-right:15px; position:relative; }

#pagetitle { padding-left:15px; background-color:#eeeeee;}

#content { background-color:#ffffff; padding-left:15px; padding-right:15px; }

/*     Styling of the round edged boxes       */
.borderleft { left:15px; position:absolute; top:0px; width:15px; height:100%; }
.borderright { right:15px; position:absolute; top:0px; width:15px; height:100%; }
.bordercentre { background-color:#ffffff; position:absolute; left:30px; right:30px; bottom:0; top:0;}
.topleftcornerborder { width:15px; height:15px; position:absolute; top:0; left:0; background-position: 0px 0px; }
.bottomleftcornerborder { width:15px; height:15px; position:absolute; bottom:0; left:0; background-position: 0px -15px; }
.toprightcornerborder { width:15px; height:15px; position:absolute; top:0; right:0; background-position: -15px 0px; }
.bottomrightcornerborder { width:15px; height:15px; position:absolute; bottom:0; right:0; background-position: -15px -15px; }
.verticalborder { background: #ffffff; position:absolute; top:15px; bottom:15px;left:0px; right:0px; }
.horizontalborder { background: #ffffff; height:15px; margin-left:15px; margin-right:15px; }