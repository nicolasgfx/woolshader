# knit_check, the command line: a program, no library. The strand demo
# (pieces_cli.cpp, --strand-demo) writes a one-strand package with no scene, so it builds the export's
# neutral scene itself: the one program with export's private headers (src/export/ARCHITECTURE.md).
knit_program(knit_check FOLDER cli
    SOURCES knit_check.cpp pieces_cli.cpp scene_cli.cpp
    DEPS model knitting fibres io settings scene jobs
    PRIVATE_ACCESS export
)

# stitchmesh_lab: the stockinette stitch-mesh experiment's command line: the
# reference surface, its analysis, the typed stitch meshes, Method A on marked pieces of an imported file (it reads
# OBJ / PLY and markings files through io), and the imported-arm commands of stitchmesh_arm.cpp (the general import,
# proposals, the reference-bound records), the yarn commands of stitchmesh_yarn.cpp and the yarn-settle commands of
# stitchmesh_settle.cpp (the knitting solver on the device: gpu). Writes OBJ, CSV, JSON and text only.
knit_program(stitchmesh_lab FOLDER cli
    SOURCES stitchmesh_lab.cpp stitchmesh_arm.cpp stitchmesh_smooth.cpp stitchmesh_yarn.cpp stitchmesh_settle.cpp
    DEPS model io knitting gpu stitchmesh jobs import
)
