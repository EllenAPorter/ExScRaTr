#include "xrtrace.h"

/**
 * Step function definition for "trace_voxel"
 */
void xrtrace_trace_voxel(cncTag_t frame, cncTag_t instance, cncTag_t i, cncTag_t j, cncTag_t k, scene_data *scene, ray_packet *rays0, ray_packet *rays1, ray_packet *rays2, ray_packet *rays3, ray_packet *rays4, ray_packet *rays5, xrtraceCtx *ctx) {

    //
    // OUTPUTS
    //

    if (i>0) {
        // Put "rays6" items
        ray_packet *rays6 = cncItemAlloc(sizeof(*rays6));
        /* TODO: Initialize rays6 */
        cncPut_rays(rays6, frame, instance+1, 0, i-1, j, k, ctx);
    }

    if (i<ctx->voxels_i-1) {
        // Put "rays7" items
        ray_packet *rays7 = cncItemAlloc(sizeof(*rays7));
        /* TODO: Initialize rays7 */
        cncPut_rays(rays7, frame, instance+1, 1, i+1, j, k, ctx);
    }

    if (j>0) {
        // Put "rays8" items
        ray_packet *rays8 = cncItemAlloc(sizeof(*rays8));
        /* TODO: Initialize rays8 */
        cncPut_rays(rays8, frame, instance+1, 2, i, j-1, k, ctx);
    }

    if (j<ctx->voxels_j-1) {
        // Put "rays9" items
        ray_packet *rays9 = cncItemAlloc(sizeof(*rays9));
        /* TODO: Initialize rays9 */
        cncPut_rays(rays9, frame, instance+1, 3, i, j+1, k, ctx);
    }

    if (k>0) {
        // Put "rays10" items
        ray_packet *rays10 = cncItemAlloc(sizeof(*rays10));
        /* TODO: Initialize rays10 */
        cncPut_rays(rays10, frame, instance+1, 4, i, j, k-1, ctx);
    }

    if (k<ctx->voxels_k-1) {
        // Put "rays11" items
        ray_packet *rays11 = cncItemAlloc(sizeof(*rays11));
        /* TODO: Initialize rays11 */
        cncPut_rays(rays11, frame, instance+1, 5, i, j, k+1, ctx);
    }

}
