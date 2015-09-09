#include "xrtrace.h"


void xrtrace_cncInitialize(xrtraceArgs *args, xrtraceCtx *ctx) {


    { // Prescribe "decompose_domain" steps
        s64 _frame;
        for (_frame = 0; _frame < ctx->num_frames; _frame++) {
            cncPrescribe_decompose_domain(_frame, ctx);
        }
    }

    { // Prescribe "camera" steps
        s64 _frame;
        for (_frame = 0; _frame < ctx->num_frames; _frame++) {
            cncPrescribe_camera(_frame, ctx);
        }
    }

    { // Prescribe "trace_voxel" steps
        s64 _frame, _instance, _i, _j, _k;
        for (_frame = 0; _frame < ctx->num_frames; _frame++) {
            for (_instance = 0; _instance < ctx->boundary_exchanges; _instance++) {
                for (_i = 0; _i < ctx->voxels_i; _i++) {
                    for (_j = 0; _j < ctx->voxels_j; _j++) {
                        for (_k = 0; _k < ctx->voxels_k; _k++) {
                            cncPrescribe_trace_voxel(_frame, _instance, _i, _j, _k, ctx);
                        }
                    }
                }
            }
        }
    }

    // Set finalizer function's tag
    xrtrace_await(ctx);

}


void xrtrace_cncFinalize(xrtraceCtx *ctx) {

}


