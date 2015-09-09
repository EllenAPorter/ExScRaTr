#include "xrtrace.h"

int cncMain(int argc, char *argv[]) {

    // Create a new graph context
    xrtraceCtx *context = xrtrace_create();

    // TODO: Set up arguments for new graph initialization
    // Note that you should define the members of
    // this struct by editing xrtrace_defs.h.
    xrtraceArgs *args = NULL;

    context->voxels_i = 2;
    context->voxels_j = 2;
    context->voxels_k = 2;
    // frames to trace
    context->num_frames = 1;
    // how many times do we expect to communicate rays
    context->boundary_exchanges = 3;

    // Launch the graph for execution
    xrtrace_launch(args, context);

    // Exit when the graph execution completes
    CNC_SHUTDOWN_ON_FINISH(context);

    return 0;
}
