#include "xrtrace.h"

/**
 * Step function defintion for "produce_image"
 */
void xrtrace_produce_image(cncTag_t frame, ray_packet *****rays, xrtraceCtx *ctx) {

    //
    // INPUTS
    //

    { // Access "rays" inputs
        s64 _neighbor, _i, _j, _k;
        for (_neighbor = 0; _neighbor < 6; _neighbor++) {
            for (_i = 0; _i < ctx->voxels_i; _i++) {
                for (_j = 0; _j < ctx->voxels_j; _j++) {
                    for (_k = 0; _k < ctx->voxels_k; _k++) {
                        /* TODO: Do something with rays[_neighbor][_i][_j][_k] */
                    }
                }
            }
        }
    }


    //
    // OUTPUTS
    //


}
