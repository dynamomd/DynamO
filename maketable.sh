#!/bin/bash
echo '<html>
<head>
<meta content="text/html;charset=ISO-8859-1" http-equiv="Content-Type">
<title></title>
</head>
<body>'

echo '<table style="text-align: left; width: 100%;" border="1" cellpadding="2" cellspacing="2">'
cat plugin.desc | gawk -F',' '{print "<tr>\n<td>\n"$1"\n</td>\n<td>\n"$2"\n</td>\n<td>\n"$3"\n</td>\n<td>"$4"\n</td>\n</tr>"}'
echo '</tbody>
</table>
</body>
</html>'