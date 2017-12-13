% MESHGEN
%
% Files
%   QuadQual       - Quality measure for quadrilateral elements.
%   TriQual        - Quality measure for triangular elements.
%   add_Edge2Elem  - Connects edges to elements.
%   add_Edges      - Add edge lists.
%   add_MLevel     - Add multilevel information.
%   add_Patches    - Adds additional patch information to the mesh.
%   dist_circ      - Signed distance for a circle.
%   dist_diff      - Signed distance for the difference between two domains.
%   dist_isect     - Signed distance for the intersection of two domains.
%   dist_rect      - Signed distance for a rectangle.
%   dist_union     - Signed distance for the union of two domains.
%   get_BdEdges    - Extract boundary edges of the mesh.
%   h_uniform      - Element size function.
%   init_Mesh      - Mesh generator.

%   load_Mesh      - Load mesh from file.
%   morph          - Morph triangular into quadrilateral meshes.
%   TProd_Mesh     - Create tensor product mesh from two one-dimensional meshes.

%   refine_REG     - Regular refinement.
%   rotate         - Rotates coordinates.
%   save_Mesh      - Save mesh to file.
%   shift          - Shifts coordinates.
%   smooth         - Laplacian smoothing.
%   stretch        - Stretches coordinates.

%   jiggle         - Mesh Jiggling.
%   dist_tri       - Signed distance for a triangle.
%   get_MeshWidth  - Computes the mesh width.

%   add_ParBd      - Add parabolic boundaries to the mesh.
%   affine_map     - generates the mapping from the reference element 
%   get_MeshMin    - Computes the mesh width.
%   graded_RefElem - generates the point set of the graded mesh in reference 
%   init_LEB       - Initialize largest edge bisection.
%   merge_Mesh     - merge Mesh2 to Mesh1 and generate a new mesh
%   orient_Elems   - check and correct the orientation of each element 
%   refine_LEB     - largest edge bisection.
%   shape_reg      - SHAP_REG shape regularity
%   tell_Mark      - lable each element from the CornerNodes set
