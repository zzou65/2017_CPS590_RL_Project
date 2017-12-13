#!/bin/csh
#
#  Purpose:
#
#    Create a GZIP'ed TAR file of the m_src/ncc_triangle files.
#
#  Modified:
#
#    02 January 2006
#
#  Author:
#
#    John Burkardt
#
#  Move to the directory just above the "ncc_triangle" directory.
#
cd $HOME/public_html/m_src
#
#  Delete any TAR or GZ file in the directory.
#
echo "Remove TAR and GZ files."
rm ncc_triangle/*.tar
rm ncc_triangle/*.gz
#
#  Create a TAR file of the "ncc_triangle" directory.
#
echo "Create TAR file."
tar cvf ncc_triangle_m_src.tar ncc_triangle/*
#
#  Compress the file.
#
echo "Compress the TAR file."
gzip ncc_triangle_m_src.tar
#
#  Move the compressed file into the "ncc_triangle" directory.
#
echo "Move the compressed file into the directory."
mv ncc_triangle_m_src.tar.gz ncc_triangle
#
#  Say goodnight.
#
echo "The ncc_triangle_m_src gzip file has been created."
