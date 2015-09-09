#include "xrtrace.h"

/**
 * Step function definition for "decompose_domain"
 */
void xrtrace_decompose_domain(cncTag_t frame, xrtraceCtx *ctx) {

    //
    // OUTPUTS
    //

    { // Put "scene" items
        s64 _i, _j, _k;
        for (_i = 0; _i < ctx->voxels_i; _i++) {
            for (_j = 0; _j < ctx->voxels_j; _j++) {
                for (_k = 0; _k < ctx->voxels_k; _k++) {
                    scene_data *scene = cncItemAlloc(sizeof(*scene));
                    /* TODO: Initialize scene */
                    cncPut_scene(scene, frame, _i, _j, _k, ctx);
                }
            }
        }
    }

}
