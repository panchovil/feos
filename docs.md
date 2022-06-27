---
project: Multicomponent Equations of State
summary: RKPR EoS
project_github: https://github.com/fedebenelli/multicomponent-eos
author: Federico Benelli
author_description: Doctoral student with focus on reservoir PVT simulation.
author_email: federico.benelli@mi.unc.edu.ar
github: https://github.com/fedebenelli
src_dir: ./src
         ./app
output_dir: ./docs
page_dir: ./docs
media_dir: ./docs/media
exclude_dir: ./example_packages
             ./test
display: public
         protected
source: true
proc_internals: true
sort: permission-alpha
print_creation_date: true
extra_mods: json_module: http://jacobwilliams.github.io/json-fortran/
            ftools: https://github.com/fedebenelli/ftools
creation_date: %Y-%m-%d %H:%M %z
md_extensions: markdown.extensions.toc
               markdown.extensions.smarty
graph: true
license: MIT
---

Fortran program for the calculation of the residual Helmholtz Energy using the
RKPR EoS (hence, PR and SRK are also included), for multicomponent mixtures.
