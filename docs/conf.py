# -- bruges documentation build configuration file --------------------------------------------------------------
#
import os
import sys
sys.path.insert(0, os.path.abspath('..'))


# Configuration file for the Sphinx documentation builder.

# -- Setup function ----------------------------------------------------------

# Defines custom steps in the process.

def autodoc_skip_member(app, what, name, obj, skip, options):
    """Exclude all private attributes, methods, and dunder methods from Sphinx."""
    import re
    exclude = re.findall(r'\._.*', str(obj))
    return skip or exclude


def remove_module_docstring(app, what, name, obj, options, lines):
    if what == "module":
        del lines[:]
    return


def setup(app):
    app.connect('autodoc-skip-member', autodoc_skip_member)
    app.connect("autodoc-process-docstring", remove_module_docstring)
    return


# -- Project information -----------------------------------------------------

project = 'bruges'
copyright = '2022, The Bruges Authors'
author = 'The Bruges Authors'


# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = 'en'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.githubpages',
              'sphinxcontrib.apidoc',
              'sphinx.ext.napoleon',
              'myst_nb',
              'sphinx.ext.viewcode',
              'matplotlib.sphinxext.plot_directive',
              'sphinx.ext.mathjax',
]

myst_enable_extensions = ["dollarmath", "amsmath"]

# MyST notebook execution
jupyter_execute_notebooks = 'force'
execution_timeout = 120

# Apidoc automation
# https://pypi.org/project/sphinxcontrib-apidoc/
# The apidoc extension and this code automatically update apidoc.
apidoc_module_dir = '../bruges'
apidoc_output_dir = './api'
apidoc_excluded_paths = []
apidoc_toc_file = False
apidoc_separate_modules = False

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']


# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_build', '_userguide']

# Entire logo SVG
# <svg width="800" height="600" version="1.1" viewBox="0 0 750 562.5" xmlns="http://www.w3.org/2000/svg">
# <g transform="translate(0 -252.4)">
# <g transform="matrix(2.635 -.8562 .8562 2.635 -21.96 -2091)" style="color-rendering:auto;color:#000000;fill:#14ca29;image-rendering:auto;isolation:auto;mix-blend-mode:normal;shape-rendering:auto;solid-color:#000000;stroke-width:1px;text-decoration-color:#000000;text-decoration-line:none;text-decoration-style:solid;text-indent:0;text-transform:none;white-space:normal">
# <path d="m-46.98 936.1q3.515-1.622 6.705-3.244 3.19-1.622 5.083-2.379 1.893-0.757 3.028-0.757 1.893 0 3.244 1.298 1.406 1.244 1.406 3.136 0 1.081-0.7029 2.271-0.6489 1.136-1.406 1.46-6.975 2.758-15.36 4.001 1.514 1.406 3.731 3.731 2.217 2.325 2.325 2.487 0.8111 1.136 2.271 2.812 1.46 1.676 2.001 2.65 0.5948 0.9192 0.5948 2.271 0 1.73-1.298 3.028-1.298 1.298-3.352 1.298-2.055 0-4.65-3.19-2.541-3.19-6.597-11.46-4.109 7.462-5.515 9.841t-2.704 3.623q-1.298 1.19-2.974 1.19-2.001 0-3.352-1.352-1.298-1.406-1.298-2.974 0-1.46 0.5407-2.217 4.975-6.759 10.38-11.73-4.542-0.7029-8.111-1.568t-7.57-2.541q-0.6489-0.3244-1.298-1.46-0.5948-1.19-0.5948-2.163 0-1.893 1.352-3.136 1.406-1.298 3.19-1.298 1.298 0 3.244 0.8111 1.947 0.757 4.921 2.271 3.028 1.46 6.867 3.298-0.7029-3.407-1.19-7.786-0.4326-4.434-0.4326-6.056 0-2.001 1.244-3.407 1.298-1.46 3.298-1.46 1.947 0 3.19 1.46 1.298 1.406 1.298 3.785 0 0.6489-0.2163 2.595-0.1622 1.893-0.5407 4.65-0.3244 2.704-0.757 6.218z" style="fill:#14ca29"/>
# </g>
# <path d="m114 275a11.08 11.08 0 0 0-10.42 9.198c-13.22 77.19-21.5 154.3-28.02 231.2-18.8 1.799-37.57 3.333-56.29 4.428a11.08 11.08 0 1 0 1.295 22.13c17.65-1.033 35.35-2.472 53.07-4.122-2.008 25.03-3.885 50.04-5.653 75.04a11.08 11.08 0 1 0 22.11 1.564c1.859-26.29 3.813-52.55 5.941-78.81 22.16-2.384 44.36-5.226 66.62-8.573 3.076 18.39 5.929 36.81 8.321 55.33a11.08 11.08 0 1 0 21.98-2.84c-2.421-18.74-5.305-37.38-8.415-55.98 18.45-3.032 36.94-6.4 55.47-10.09-4.171 3.89-7.792 7.763-10.87 11.59-12.06 14.99-17.45 30.09-11.09 43.2 3.183 6.55 9.938 11.21 17.01 12.44 7.068 1.226 14.39-0.1838 21.94-3.436 10.44-4.496 21.58-12.76 33.11-25.35 11.22 78.97-11.75 150.1-41.57 188.1-16.38 20.83-34.06 30.54-47.05 29.66-12.98-0.8812-26.23-11.03-36.1-41.82a11.08 11.08 0 1 0-21.11 6.765c11.34 35.38 31.62 55.53 55.71 57.17 24.08 1.634 46.99-13.93 65.97-38.07 37.96-48.28 63.06-135.4 41.56-230.1a11.08 11.08 0 0 0-19.85-3.95c-15.5 21.89-29.99 33.21-39.44 37.28-4.726 2.036-8.047 2.186-9.384 1.954s-0.7451-0.0477-0.8609-0.2857c-0.2316-0.4765-0.7676-8.197 8.419-19.61 9.186-11.42 26.68-25.66 55.05-39.08a11.08 11.08 0 0 0-7.289-20.8c-37.95 8.973-75.73 16.47-113.4 22.64-13.48-74.08-32.28-147.1-55.78-219.1a11.08 11.08 0 0 0-11.03-7.63zm4.009 59.86c16.67 56.19 30.44 112.9 40.9 170.3-20.34 3.053-40.63 5.645-60.88 7.881 5.091-59.48 11.24-118.9 19.98-178.2zm326.7 3.282a11.08 11.08 0 0 0-10.94 7.293c-10.22 26.98-33.58 91.48-34.64 152.6-0.5292 30.58 4.485 60.8 21.13 84.77 16.64 23.97 45.1 40.28 86.12 42.97 15.68 1.028 29.67 0.4913 42.2-1.223 34.18 50.22 102 89.36 150.3 111.8a11.08 11.08 0 1 0 9.342-20.1c-43.72-20.32-104.6-56.93-136.1-97.03 13-4.309 23.53-10.37 31.16-18.07 8.908-9 13.76-20.83 12.41-32.65-1.353-11.82-8.486-22.64-19.47-31.12-13.27-10.25-25.74-15.33-37.71-14.09-11.96 1.232-21.44 9.869-25.83 20.17-7.142 16.73-5.358 38.48 4.07 60.03-8.674 0.7494-18.33 0.8804-28.97 0.1829-36.1-2.367-56.58-15.09-69.36-33.5-12.78-18.41-17.65-43.91-17.17-71.74 0.9633-55.67 23.22-118.8 33.21-145.2a11.08 11.08 0 0 0-9.788-15.14zm-93.08 105.6a11.08 11.08 0 0 0-7.721 19.03l1.173 1.173a11.08 11.08 0 1 0 15.67-15.67l-1.173-1.173a11.08 11.08 0 0 0-7.952-3.358zm3.948 49.46a11.08 11.08 0 0 0-10.85 11.3l0.6952 126.6a11.08 11.08 0 1 0 22.16-0.1219l-0.6952-126.6a11.08 11.08 0 0 0-11.31-11.18zm205.2 52.15c3.713-0.3825 11.26 1.375 21.89 9.586 7.549 5.832 10.47 11.47 11 16.1s-0.9047 9.245-6.143 14.54c-5.258 5.312-14.67 10.55-28.08 14.1-9.317-19.06-10.19-38.58-6.379-47.5 2.142-5.018 4.007-6.442 7.72-6.824z" style="color-rendering:auto;color:#000000;fill-rule:evenodd;fill:#7f7f7f;image-rendering:auto;isolation:auto;mix-blend-mode:normal;shape-rendering:auto;solid-color:#000000;text-decoration-color:#000000;text-decoration-line:none;text-decoration-style:solid;text-indent:0;text-transform:none;white-space:normal"/>
# </g>
# </svg>

html_theme_options = {
    "sidebar_hide_name": True,
    "footer_icons": [
        {
            "name": "Agile",
            "url": "https://code.agilescientific.com",
            "html": """
                <svg width="200" height="200" version="1.1" viewBox="0 0 187.5 187.5" xmlns="http://www.w3.org/2000/svg">
                <g transform="matrix(3.94 -1.28 1.28 3.94 -906.1 -3676)" style="color-rendering:auto;color:#000000;fill:#14ca29;image-rendering:auto;isolation:auto;mix-blend-mode:normal;shape-rendering:auto;solid-color:#000000;stroke-width:1px;text-decoration-color:#000000;text-decoration-line:none;text-decoration-style:solid;text-indent:0;text-transform:none;white-space:normal">
                <path d="m-46.98 936.1q3.515-1.622 6.705-3.244t5.083-2.379 3.028-0.757q1.893 0 3.244 1.298 1.406 1.244 1.406 3.136 0 1.081-0.7029 2.271-0.6489 1.136-1.406 1.46-6.975 2.758-15.36 4.001 1.514 1.406 3.731 3.731t2.325 2.487q0.8111 1.136 2.271 2.812 1.46 1.676 2.001 2.65 0.5948 0.9192 0.5948 2.271 0 1.73-1.298 3.028t-3.352 1.298q-2.055 0-4.65-3.19-2.541-3.19-6.597-11.46-4.109 7.462-5.515 9.841t-2.704 3.623q-1.298 1.19-2.974 1.19-2.001 0-3.352-1.352-1.298-1.406-1.298-2.974 0-1.46 0.5407-2.217 4.975-6.759 10.38-11.73-4.542-0.7029-8.111-1.568t-7.57-2.541q-0.6489-0.3244-1.298-1.46-0.5948-1.19-0.5948-2.163 0-1.893 1.352-3.136 1.406-1.298 3.19-1.298 1.298 0 3.244 0.8111 1.947 0.757 4.921 2.271 3.028 1.46 6.867 3.298-0.7029-3.407-1.19-7.786-0.4326-4.434-0.4326-6.056 0-2.001 1.244-3.407 1.298-1.46 3.298-1.46 1.947 0 3.19 1.46 1.298 1.406 1.298 3.785 0 0.6489-0.2163 2.595-0.1622 1.893-0.5407 4.65-0.3244 2.704-0.757 6.218z" style="fill:#14ca29"/>
                </g>
                </svg>
            """,
            "class": "",
        },
    ],
}

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'furo'
html_logo = '_static/bruges_logo.png'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
html_favicon = '_static/favicon.ico'


html_css_files = [
    'custom.css',
]
