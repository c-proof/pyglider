import yaml
import glob
import os
import xarray as xr


"""
Utilities to make smart directories for the websites.  Dumb persons erdapp
server...
"""

header = """
<!DOCTYPE html>
<html lang="en">

<head>

<meta charset="utf-8">
<meta http-equiv="X-UA-Compatible" content="IE=edge">
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="description" content="Canadian-Pacific Robotic Ocean Observing Facility">
<meta property="og:title" content="Data" />
<meta property="og:type" content="website" />
<meta property="og:url" content="http://cproof.uvic.ca///data/" />
<meta property="og:image" content="http://cproof.uvic.ca//img/ControlRoomMar19.jpg">



<link rel='shortcut icon' type='image/png' href='/favicon.png' />

<title>Data - C-PROOF</title>

<link rel="canonical" href="http://cproof.uvic.ca/data/">

<!-- Bootstrap Core CSS -->
<link rel="stylesheet" href="http://cproof.uvic.ca/css/bootstrap.css">

<!-- Custom CSS -->
<link rel="stylesheet" href="http://cproof.uvic.ca/css/clean-blog.css">

<!-- Adjust Colors -->
<link rel="stylesheet" href="http://cproof.uvic.ca/colorscheme.css">

<!-- Pygments Github CSS -->
<link rel="stylesheet" href="http://cproof.uvic.ca/css/syntax.css">

<!-- Custom Fonts -->
<link href="https://maxcdn.bootstrapcdn.com/font-awesome/4.2.0/css/font-awesome.min.css" rel="stylesheet" type="text/css">
<link href='https://fonts.googleapis.com/css?family=Lora:400,700,400italic,700italic' rel='stylesheet' type='text/css'>
<link href='https://fonts.googleapis.com/css?family=Open+Sans:300italic,400italic,600italic,700italic,800italic,400,300,600,700,800' rel='stylesheet' type='text/css'>

<!-- HTML5 Shim and Respond.js IE8 support of HTML5 elements and media queries -->
<!-- WARNING: Respond.js doesn't work if you view the page via file:// -->
<!--[if lt IE 9]>
    <script src="https://oss.maxcdn.com/libs/html5shiv/3.7.0/html5shiv.js"></script>
    <script src="https://oss.maxcdn.com/libs/respond.js/1.4.2/respond.min.js"></script>
<![endif]-->
<!-- Math Jax -->
<script type="text/x-mathjax-config">
MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}});
</script>
<script type="text/javascript"
src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
</script>
<script type="text/x-mathjax-config">
MathJax.Hub.Config({
    TeX: { equationNumbers: { autoNumber: "AMS" } }
});
</script>

<!--<script src="http://cdnjs.cloudflare.com/ajax/libs/d3/3.4.13/d3.min.js" type="text/javascript"></script>-->

<!-- jQuery -->
<script src="https://code.jquery.com/jquery-2.1.4.min.js"></script>

<!-- Bootstrap Core JavaScript -->
<script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/js/bootstrap.min.js"></script>

<!-- Custom Theme JavaScript -->
<script src="http://cproof.uvic.ca/js/clean-blog.min.js "></script>



</head>
<body>

<!-- Navigation -->
<nav class="navbar navbar-default navbar-custom navbar-fixed-top">
<div class="container-fluid">
    <!-- Brand and toggle get grouped for better mobile display -->
    <div class="navbar-header page-scroll">
        <button type="button" class="navbar-toggle" data-toggle="collapse" data-target="#bs-example-navbar-collapse-1">
            <span class="sr-only">Toggle navigation</span>
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
        </button>
        <a class="navbar-brand" href="/">C-PROOF</a>
    </div>

    <!-- Collect the nav links, forms, and other content for toggling -->
    <div class="collapse navbar-collapse" id="bs-example-navbar-collapse-1">
        <ul class="nav navbar-nav navbar-right">
            <li>
                <a href="/">Home</a>
            </li>

            <li>
                <a href="/about.html">About</a>
            </li>

            <li>
                <a href="/platforms/">Platforms</a>
            </li>

            <li>
                <a href="/data/">Data</a>
            </li>

            <li>
                <a href="/deployments/">Deployments</a>
            </li>

            <li>
                <a href="/members.html">People</a>
            </li>

        </ul>
    </div>
    <!-- /.navbar-collapse -->
</div>
<!-- /.container -->
</nav>


<!-- Page Header -->
<header class="intro-header" style="background-image: url('/img/ControlRoomMar19.jpg')">
<div class="container">
    <div class="row">
        <div class="col-lg-8 col-lg-offset-2 col-md-10 col-md-offset-1">
            <div class="site-heading" style="padding: 75px 0">

                <h1>Data</h1>
                <hr class="small">
                <span class="subheading">C-PROOF Data</span>
            </div>
        </div>
    </div>
</div>
</header>
"""


def index_deployments(dir):
    """
    Get useful info from deployments under "dir" and add to an
    index.html in "dir"

    The structure is meant to be simple:

    dir/deployment1/deployment.yml
    dir/deployment2/deployment.yml
    dir/deployment3/deployment.yml

    Makes a file `dir/index.html` that has table of the deployments.
    """
    subdirs = glob.glob(dir + '/*')
    html = header
    html += """
    <table border="1">
    <tr><th>Deployment</th><th>Description</th><th>Dates</th><th>Geographic Extent</th></tr>\n"""
    for d in subdirs:
        if os.path.isdir(d):
            if 1:
                nc = glob.glob(d+'/L1-timeseries/*.nc')
                with xr.open_dataset(nc[0]) as ds:
                    att = ds.attrs
                    print(att)
                    dname = att['deployment_name']
                    html += '    <tr>'
                    html += f'<td><a href="{dname}">{dname}</a></td>'
                    html += f"<td>{att['comment']}</td> "
                    html += f" <td>{att['time_coverage_start']}<br/>{att['time_coverage_end']}</td>"
                    html += f" <td>{att['geospatial_lat_min']:3.2f} {att['geospatial_lat_max']:3.2f} N<br/>{att['geospatial_lon_min']:3.2f} {att['geospatial_lon_max']:3.2f} W</td>"
                    html += '    </tr>'

            else:
                pass
            html += '\n'
    html += "</table></html>"
    with open(dir + './index.html', 'w') as fout:
        fout.write(html)
