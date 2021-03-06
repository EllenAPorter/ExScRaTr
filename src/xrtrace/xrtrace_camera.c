#include "xrtrace.h"

/**
 * Step function definition for "camera"
 */
void xrtrace_camera(cncTag_t frame, xrtraceCtx *ctx) {

    //
    // OUTPUTS
    //

    { // Put "rays" items
        s64 _neighbor, _i, _j, _k;
        for (_neighbor = 0; _neighbor < 6; _neighbor++) {
            for (_i = 0; _i < ctx->voxels_i; _i++) {
                for (_j = 0; _j < ctx->voxels_j; _j++) {
                    for (_k = 0; _k < ctx->voxels_k; _k++) {
                        ray_packet *rays = cncItemAlloc(sizeof(*rays));
                        /* TODO: Initialize rays */
                        cncPut_rays(rays, frame, 0, _neighbor, _i, _j, _k, ctx);
                    }
                }
            }
        }
    }

}
