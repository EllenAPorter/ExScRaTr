#include <cnc/cnc.h>
#include <cnc/debug.h>

#include "util.h"

using namespace CnC;

struct ray_tracer;

typedef std::pair<int,int> pair;

/**
 * Consumes a set of pixel locations (2D)
 * and ray traces the scene for them
 */
struct step_trace_ray {
	int execute ( const pair & frame_fragment, ray_tracer & rt ) const;
};

/**
 * Composes a single image from image fragments
 */
struct step_compose_image {
	int execute ( const int & frame, ray_tracer & rt ) const;
};

/**
 * Dynamically determines when to rearrange the pixel
 * allocations for the ray tracer step
 */
struct step_decompose_image {
	int execute ( const int & frame, ray_tracer & rt ) const;
};

/**
 * Produces a new dynamic scene for every frame
 */
struct step_compute_dynamic_scene {
	int execute ( const int & frame, ray_tracer & rt ) const;
};

/**
 * CnC Ray Tracer context
 */
struct ray_tracer : public context< ray_tracer >
{
	// Step collections
	step_collection < step_trace_ray > ray_trace;
	step_collection < step_compose_image > compose_image;
	step_collection < step_decompose_image > decompose_image;
	step_collection < step_compute_dynamic_scene > compute_dynamic_scene;

	// Item collections
	item_collection < int, std::vector<Primitive*> > dynamic_scene;
	item_collection < pair, std::vector<Point2D*> > pixel_locations;
	item_collection < pair, std::vector<Pixel*> > pixels;

	item_collection< pair, double > execution_time;

	// Tag collections
	tag_collection < int > frame;
	tag_collection < pair > frame_fragment;

	// Constants
	int number_frames;
	int number_fragments;

	int image_height;
	int image_width;

	int number_pixel_samples;
	int number_light_samples;

	std::vector<Primitive*> static_scene;
	std::vector<Luminaire*> luminaires;

	ray_tracer(): context< ray_tracer >(),

			// Initialize step collections
			ray_trace ( *this, "step ray trace" ),
			compose_image ( *this, "step compose image" ),
			decompose_image ( *this, "step decompose image" ),
			compute_dynamic_scene( *this, "step compute dynamic scene"),

			// Initialize item collections
			pixel_locations ( *this, "pixel locations" ),
			dynamic_scene ( *this, "dynamic scene objects" ),
			pixels ( *this, "pixels with their locations" ),
			execution_time( *this, "execution time" ),

			// Initialize tag collections
			frame( *this, "per frame"),
			frame_fragment ( *this, "per frame per fragments" )
			{
					// Tag control
					compose_image.controls(frame);
					decompose_image.controls(frame_fragment);

					// Step generation
					frame_fragment.prescribes(ray_trace, *this );

					frame.prescribes(compute_dynamic_scene, *this );
					frame.prescribes(decompose_image, *this );
					frame.prescribes(compose_image, *this );

					// Producers
					ray_trace.produces(pixels);
					decompose_image.produces(pixel_locations);
					compute_dynamic_scene.produces(dynamic_scene);

					// Consumers
					ray_trace.consumes(pixel_locations);
					ray_trace.consumes(dynamic_scene);

					compose_image.consumes(pixels);
			};
};
