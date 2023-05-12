### usage vmd -e visualise_pathway_hole.tcl -args  Zhang_7dwb.pdb Zhang_7dwb.vmd ###

set FILE_STRUCTURE [lindex $argv 0]
set FILE_PORE_SURFACE [lindex $argv 1]
puts $FILE_STRUCTURE
puts $FILE_PORE_SURFACE

mol new $FILE_STRUCTURE

#axes location Off
color Display Background white
display projection Orthographic
display depthcue off
display ambientocclusion on
display rendermode GLSL

mol delrep 0 top
# protein:
mol addrep top
mol modselect 0 top (protein)
mol modstyle 0 top NewRibbons 0.3 10.0 4.1 0
mol modcolor 0 top ColorID 6
mol modmaterial 0 top AOEdgy

### Load surface ###
echo $FILE_PORE_SURFACE
source $FILE_PORE_SURFACE
#source "Zhang_7dwb.vmd"
