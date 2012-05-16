<?php header("Content-type: text/css"); 
$sitebg = "#b7c7c8";
$pagebg = "#ffffff";
$leftmenuhovercolor = "#eeeeee";
$codebg = "#ccdddd";
$inactivebutton = "#ccdddd";
$hoverbutton = "#bbcccc";
?>

/* Main html element styles */
html { height: 100%; }
/*#e4f4ff*/
body { 
    margin: 0; 
    padding: 0; 
    height: 100%; 
    background-color:<?php echo $sitebg; ?>; 
}

/* Floating footer code */
#wrapper { 
    min-height:100%; 
    min-width:600px;
    overflow-y:visible; 
    position:relative; 
}

#wrapperfooterpad { 
    height:55px; 
    margin-bottom:10px; 
    clear:both; 
    position:relative;
    margin-left:190px; 
}

#header { 
    position:relative; 
    padding-left:30px; 
    padding-right:30px; 
    height:145px;
}

#sitelogo { 
    background-repeat:no-repeat; 
    background-image:url(../images/sitelogo.png); 
    background-position: 15px 15px; 
}

#headercontentspacing { 
    height: 15px; 
    clear: both; 
}

#footer {
    height:40px; 
    margin-top:-55px;
    position:relative;
}

.w3footerlogo { 
    border:0; 
    width:88px; 
    height:31px; 
    position:absolute; 
    top:5px; 
}

.sprite { 
    background-image: url(../images/csssprites.png); 
}

#leftmenu { 
    width:160px; 
    margin-left:15px; 
    position:absolute;
}

#leftmenu a { 
    color:#000000; 
    text-decoration: none; 
    line-height:30px; 
    display:block; 
    padding-left:16px; 
    background-color:<?php echo $pagebg; ?>;
}

#leftmenu a:hover { 
    background-color:<?php echo $leftmenuhovercolor;?>; 
}


#contentwrapper { 
    margin-left:190px; 
    margin-right:15px;
    position:relative;
}

#pagetitle { 
    position:absolute; 
    font-weight:bold; 
    top:2px; 
    padding-left:15px;
}

#content { 
    background-color:<?php echo $pagebg; ?>; 
    padding-left:15px; 
    padding-right:15px; 
    padding-top:10px; 
    padding-bottom:10px;
}

.code { background-color:<?php echo $codebg;?>; font-family:monospace; padding:5px; }

/*     Styling of the round edged boxes       */
.borderleft { left:15px; position:absolute; top:0px; width:15px; height:100%; }
.borderright { right:15px; position:absolute; top:0px; width:15px; height:100%; }
.bordercentre { background-color:<?php echo $pagebg; ?>; position:absolute; left:30px; right:30px; bottom:0; top:0;}
.verticalborder { background: <?php echo $pagebg; ?>; position:absolute; top:15px; bottom:15px;left:0px; right:0px; }
.horizontalborder { background: <?php echo $pagebg; ?>; height:15px; margin-left:15px; margin-right:15px; }

.topleftcornerborder { background-image: url(../images/csssprites.png); width:15px; height:15px; position:absolute; top:0; left:0; background-position: 0 0; }
.bottomleftcornerborder { background-image: url(../images/csssprites.png); width:15px; height:15px; position:absolute; bottom:0; left:0; background-position: 0 -15px; }
.toprightcornerborder { background-image: url(../images/csssprites.png); width:15px; height:15px; position:absolute; top:0; right:0; background-position: -15px 0; }
.bottomrightcornerborder { background-image: url(../images/csssprites.png); width:15px; height:15px; position:absolute; bottom:0; right:0; background-position: -15px -15px; }

/* Button */
.button { position:relative; padding-left:15px; padding-right:15px; height:30px; display:inline-block; }
.button span { position:absolute; width:100%; height:100%; top:0; left:0; z-index:1; }
.button .center { background-color: <?php echo $inactivebutton;?>; height:30px; line-height:30px; display:inline-block; }
.button .left { background-image: url(../images/csssprites.png); width:15px; height:30px; position:absolute; top:0; left:0; background-position: -35px 0; }
.button .right { background-image: url(../images/csssprites.png); width:15px; height:30px; position:absolute; bottom:0; right:0; background-position: -50px 0; }

.button:hover .center { background-color: <?php echo $hoverbutton;?>; }
.button:hover .left { background-position: -70px 0; }
.button:hover .right { background-position: -85px 0; }

.bordercentre p { padding:0; margin:0; line-height:40px; }

h1 { font-size:18pt; border-bottom:3px solid; display:block; clear:both; }

h2 { font-size:16pt; display:block; }

.figure { float: right; display: inline-block; margin:5px; border: solid 1px; }
.figure iframe { width:100%; }
.figure .caption { text-align:center; clear: both;padding:5px; }
