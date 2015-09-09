/******************************************************************************
// * CnC exascale ready ray tracer
// *
// * Author: Ellen Porter (ellen.porter@wsu.edu)
// *         Washington State University
// *
//****************************************************************************/

/******************************************************************************
// * Graph parameters */
$context {
	
	// decomposed domain size
	int voxels_i;
	int voxels_j;
	int voxels_k;
	
	// frames to trace
	int num_frames;
	
	// how many times do we expect to communicate rays
	int boundary_exchanges;
	
};

/******************************************************************************
//* Item collection declarations */

// scene data for each voxel, produced by decompose_domain
[ scene_data *scene : frame, i, j, k ];

// ray data passed to neighbors, produced by camera and trace_voxel
[ ray_packet *rays : frame, instance, neighbor, i, j, k ];

/******************************************************************************
//* CnC steps */

( decompose_domain : frame )
-> [ scene : frame, $range(0, #voxels_i), $range(0, #voxels_j), $range(0, #voxels_k) ];

( camera : frame )
-> [ rays : frame, 0, 2, $range(0, #voxels_i), 0, $range(0, #voxels_k) ];

( trace_voxel : frame, instance, i, j, k )
<- [ scene: frame, i, j, k],
   [ rays : frame, instance  , 0, i  , j  , k   ] $when(instance > 0 && i<#voxels_i-1),
   [ rays : frame, instance  , 1, i  , j  , k   ] $when(instance > 0 && i>0),
   [ rays : frame, instance  , 2, i  , j  , k   ] $when(instance == 0 && j == 0 || j<#voxels_j-1),
   [ rays : frame, instance  , 3, i  , j  , k   ] $when(instance > 0 && j>0),
   [ rays : frame, instance  , 4, i  , j  , k   ] $when(instance > 0 && k<#voxels_k-1),
   [ rays : frame, instance  , 5, i  , j  , k   ] $when(instance > 0 && k>0)
-> [ rays : frame, instance+1, 0, i-1, j  , k   ] $when(i>0),
   [ rays : frame, instance+1, 1, i+1, j  , k   ] $when(i<#voxels_i-1),
   [ rays : frame, instance+1, 2, i  , j-1, k   ] $when(j>0),
   [ rays : frame, instance+1, 3, i  , j+1, k   ] $when(j<#voxels_j-1),
   [ rays : frame, instance+1, 4, i  , j  , k-1 ] $when(k>0),
   [ rays : frame, instance+1, 5, i  , j  , k+1 ] $when(k<#voxels_k-1);


/******************************************************************************
//* Input output relationships from environment */
( $initialize: () )
-> (decompose_domain : $range(0, #num_frames)),
   (camera           : $range(0, #num_frames)),
   (trace_voxel      : $range(0, #num_frames), $range(0, #boundary_exchanges), $range(0, #voxels_i), $range(0, #voxels_j), $range(0, #voxels_k));

( $finalize: () );
