#include <X11/Xlib.h>

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>
#include <assert.h>

#include <math.h>
#include <string.h>

#include <poll.h>
#include <time.h>

#include <pthread.h>
#include <semaphore.h>

#define MAX_FPS 15

#define WIDTH  800
#define HEIGHT 800

#define WINDOW_BORDER_WIDTH 1


// ** Pixels **
typedef uint32_t Pixel;

#define PIXEL_DEPTH 24

// hard coded for Little-Endian
#define PIXEL_RGB(r, g, b) (Pixel) ( \
	(r << 8 * 2) |					 \
	(g << 8    ) |					 \
	(b		   )  					 \
)
#define PIXEL_R(p) ((p >> 8 * 2) & 0xff)
#define PIXEL_G(p) ((p >> 8    ) & 0xff)
#define PIXEL_B(p) ((p		   ) & 0xff)

// ** Pixel Buffer **
Pixel pixel_buffer[WIDTH * HEIGHT] = {0};
	


// ** Vectors **
typedef struct V3 {
	float x;
	float y;
	float z;
} V3;

float dot_product(V3 a, V3 b) {
	return (
		a.x * b.x +
	   	a.y * b.y + 
		a.z * b.z
	);
}

float norm(V3 a) {
	return sqrt(
		a.x * a.x +
		a.y * a.y +
		a.z * a.z
	);
}

V3 normalize(V3 a) {
	float n = norm(a);

	return (V3) {
		a.x / n,
		a.y / n,
		a.z / n
	};
}

V3 scalar_product(V3 a, float b) {
	return (V3) {
		b * a.x,
		b * a.y,
		b * a.z
	};
}

V3 vector_sum(V3 a, V3 b) {
	return (V3) {
		a.x + b.x,
		a.y + b.y,
		a.z + b.z
	};
}
V3 vector_sub(V3 a, V3 b) {
	return (V3) {
		a.x - b.x,
		a.y - b.y,
		a.z - b.z
	};
}

#define VECTOR_CMP_ERROR_TOLERANCE 0.00001
int vector_cmp(V3 a, V3 b) {
	return (
		fabs(a.x - b.x) <= VECTOR_CMP_ERROR_TOLERANCE &&
		fabs(a.y - b.y) <= VECTOR_CMP_ERROR_TOLERANCE &&
		fabs(a.z - b.z) <= VECTOR_CMP_ERROR_TOLERANCE
	);
}

// check if a is between b and c
int point_in_segment(V3 a, V3 b, V3 c) {
	V3 ba = vector_sub(a, b);
	V3 bc = vector_sub(c, b);

	return (
		(vector_cmp(normalize(ba), normalize(bc))) &&
		(norm(ba) <= norm(bc))
	);
}

// ** World **
typedef struct Sphere {
	float radius;
} Sphere;

typedef struct Plane {
	V3 normal;
} Plane;

typedef struct LightSource {
	Pixel colour;
	float intensity;
} LightSource;

typedef struct Camera {
	// the plane has its center on the local x-axis 
	// and it's orhogonal to it
	float near_plane_distance; // relative to eye
	float near_plane_width;
	float near_plane_height;

	float far_plane_distance;
} Camera;


typedef enum ObjectType {
	// primitive shapes
	SPHERE,
	PLANE,

	// special objects
	LIGHT_SOURCE,
	CAMERA,
} ObjectType;

typedef struct Object {
	ObjectType type;
	V3 position;

	union {
		Sphere sphere;
		Plane plane;
		LightSource light_source;
		Camera camera;
	};
} Object;


// ** Ray tracer **
typedef struct Ray {
	V3 origin;
	V3 direction;
} Ray;

#define MAX_OBJECT_INTERSECTIONS 2
int ray_object_intersect(Object object, Ray ray, V3* intersections) {
	assert(intersections);

	// normalize ray direction vector
	ray.direction = normalize(ray.direction);

	switch(object.type) {
		case SPHERE: 
			V3 _oc = vector_sub(ray.origin, object.position);

			float _alpha = dot_product(ray.direction, _oc);
			float _beta  = (dot_product(_oc, _oc) - object.sphere.radius * object.sphere.radius);
			float _delta = _alpha * _alpha - _beta;

			if(_delta < 0) {
				return 0; // there are no intersections
			}

			float t1 = -_alpha + sqrt(_delta);
			float t2 = -_alpha - sqrt(_delta);

			int ret = 0;
			if(t1 >= 0) {
				intersections[ret++] = vector_sum(ray.origin, scalar_product(ray.direction, t1));
			}
			if(t2 >= 0) {
				intersections[ret++] = vector_sum(ray.origin, scalar_product(ray.direction, t2));
			}

			return ret;
		case PLANE: 
			float cos_au = dot_product(object.plane.normal, ray.direction);
			if(!cos_au) {
				// plane and ray are parallel
				return 0;
			}

			float cos_aoo = dot_product(object.plane.normal, vector_sub(object.position, ray.origin));
			float t = cos_aoo / cos_au;

			if(t < 0) {
				return 0;
			}

			intersections[0] = vector_sum(ray.origin, scalar_product(ray.direction, t));
			return 1;
		default:
	}

	return -1;
}

int render_scene(Object* scene, size_t scene_size, Pixel* out_buffer, size_t width, size_t height) {
	assert(scene);
	assert(out_buffer);
	assert(width  > 0);
	assert(height > 0);

	// search for camera and light source in the scene
	LightSource light;
	V3 light_position;
	Camera camera;
	V3 camera_position;

	int found_light, found_camera;
	for(size_t i = 0; i < scene_size; i++) {
		if(scene[i].type == LIGHT_SOURCE) {
			light = scene[i].light_source;
			light_position = scene[i].position;
			found_light = 1;
		} else 
		if(scene[i].type == CAMERA) {
			camera = scene[i].camera;
			camera_position = scene[i].position;
			found_camera = 1;
		}
	}
	assert(found_light && found_camera);

	// TODO: world coordinates
	// for now we just assert that the camera position is (0, 0, 0), 
	// because we haven't figured out actual world coordinates yet.
	assert(!(
		camera_position.x ||
		camera_position.y ||
		camera_position.z
	));


	// boundary planes
	V3 near_plane_center = {
		camera.near_plane_distance,
		0,
		0
	};

	// for every pixel (i, j)
	for(size_t i = 0; i < height; i++) {
		for(size_t j = 0; j < width; j++) {
			out_buffer[i * width + j] = PIXEL_RGB(0x00, 0x00, 0x00);

			// figure out where (i, j) is on the near plane
			float near_plane_i = (float)i / height * camera.near_plane_height;
			float near_plane_j = (float)j / width  * camera.near_plane_width;
			float near_plane_x = near_plane_j - camera.near_plane_width  / 2.0;
			float near_plane_y = near_plane_i - camera.near_plane_height / 2.0;

			// ray 
			Ray ray = {
				.origin = camera_position,
				.direction = (V3){
					near_plane_center.x,
					near_plane_x,
					near_plane_y
				},
			};

			// find closest intersection 
			ssize_t closest_intersection_object = -1;
			V3 closest_intersection = {
				camera.far_plane_distance,
				0,
				0,
			};
			for(size_t obj = 0; obj < scene_size; obj++) {
				if(scene[obj].type == SPHERE || scene[obj].type == PLANE) {
					V3 intersections[MAX_OBJECT_INTERSECTIONS];
					int n_intersections = ray_object_intersect(scene[obj], ray, intersections);
					assert(n_intersections >= 0);
					assert(n_intersections <= MAX_OBJECT_INTERSECTIONS);

					for(size_t intersection = 0 ; intersection < (size_t)n_intersections; intersection++) {
						if(closest_intersection.x > intersections[intersection].x) {
							closest_intersection_object = obj;
							closest_intersection = intersections[intersection];
						}
					}
				}
			}
			if(closest_intersection_object < 0) continue;

			V3 intersection_normal;
			switch(scene[closest_intersection_object].type) {
				case SPHERE: 
					intersection_normal = normalize(vector_sub(closest_intersection, scene[closest_intersection_object].position)); 
					break;
				case PLANE:  
					intersection_normal = scene[closest_intersection_object].plane.normal; 
					break;
				default: 
					assert(0 && "Unreachable");
			}

			if(dot_product(intersection_normal, normalize(ray.direction)) > 0.0) {
				// the surface is facing away
				continue;
			}

			// find intersection -> light source vector
			Ray ray1 = {
				.direction = vector_sub(light_position, closest_intersection),
				.origin = closest_intersection,
			};

			// check for occlusions (shadows)
			int occluded = 0;
			for(size_t obj = 0; obj < scene_size; obj++) {
				if(obj != (size_t)closest_intersection_object && (scene[obj].type == SPHERE || scene[obj].type == PLANE)) {
					V3 intersections[MAX_OBJECT_INTERSECTIONS];
					int n_intersections = ray_object_intersect(scene[obj], ray1, intersections);
					assert(n_intersections >= 0);
					assert(n_intersections <= MAX_OBJECT_INTERSECTIONS);

					for(size_t intersection = 0 ; intersection < (size_t)n_intersections; intersection++) {
						if(point_in_segment(intersections[intersection], closest_intersection, light_position)) {
							// it is an occlusion
							occluded = 1;
							break;
						}
					}
				}
				if(occluded) break;
			}
			if(occluded) continue;
		
			float reflectance = dot_product(intersection_normal, normalize(ray1.direction));
			// assert(reflectance <= 1.0);
			
			if(reflectance < 0.0) {
				// the surface is facing away from the light source
				continue;
			}

			// assert(light.intensity <= 1.0 && light.intensity >= 0.0);
			Pixel p = PIXEL_RGB(
				(Pixel) (PIXEL_R(light.colour) * reflectance * light.intensity),
				(Pixel) (PIXEL_G(light.colour) * reflectance * light.intensity),
				(Pixel) (PIXEL_B(light.colour) * reflectance * light.intensity)
			);

			out_buffer[i * width + j] = p;
		}
	}
	

	return 0;
}


int render_scene_threaded(Object* scene, size_t scene_size, Pixel* out_buffer, size_t width, size_t height, int n_threads, int thread) {
	assert(scene);
	assert(out_buffer);
	assert(width  > 0);
	assert(height > 0);

	// search for camera and light source in the scene
	LightSource light;
	V3 light_position;
	Camera camera;
	V3 camera_position;

	int found_light, found_camera;
	for(size_t i = 0; i < scene_size; i++) {
		if(scene[i].type == LIGHT_SOURCE) {
			light = scene[i].light_source;
			light_position = scene[i].position;
			found_light = 1;
		} else 
		if(scene[i].type == CAMERA) {
			camera = scene[i].camera;
			camera_position = scene[i].position;
			found_camera = 1;
		}
	}
	assert(found_light && found_camera);

	// TODO: world coordinates
	// for now we just assert that the camera position is (0, 0, 0), 
	// because we haven't figured out actual world coordinates yet.
	assert(!(
		camera_position.x ||
		camera_position.y ||
		camera_position.z
	));


	// boundary planes
	V3 near_plane_center = {
		camera.near_plane_distance,
		0,
		0
	};


	// TODO: support more thread configurations
	assert(n_threads == 16); 

	assert(thread >= 0);
	assert(thread < n_threads);

	size_t grid_w = 4;
	size_t grid_h = 4;

	// assuming that you can divide the screen by 4
	size_t grid_cell_w = width  / grid_w;
	size_t grid_cell_h = height / grid_h;
	
	size_t i0 = (thread / grid_w) * grid_cell_h;
	size_t j0 = (thread % grid_w) * grid_cell_w;

	// for every pixel (i, j)
	for(size_t i = i0; i < i0 + grid_cell_h; i++) {
		for(size_t j = j0; j < j0 + grid_cell_w; j++) {
			out_buffer[i * width + j] = PIXEL_RGB(0x00, 0x00, 0x00);

			// figure out where (i, j) is on the near plane
			float near_plane_i = (float)i / height * camera.near_plane_height;
			float near_plane_j = (float)j / width  * camera.near_plane_width;
			float near_plane_x = near_plane_j - camera.near_plane_width  / 2.0;
			float near_plane_y = near_plane_i - camera.near_plane_height / 2.0;

			// ray 
			Ray ray = {
				.origin = camera_position,
				.direction = (V3){
					near_plane_center.x,
					near_plane_x,
					near_plane_y
				},
			};

			// find closest intersection 
			ssize_t closest_intersection_object = -1;
			V3 closest_intersection = {
				camera.far_plane_distance,
				0,
				0,
			};
			for(size_t obj = 0; obj < scene_size; obj++) {
				if(scene[obj].type == SPHERE || scene[obj].type == PLANE) {
					V3 intersections[MAX_OBJECT_INTERSECTIONS];
					int n_intersections = ray_object_intersect(scene[obj], ray, intersections);
					assert(n_intersections >= 0);
					assert(n_intersections <= MAX_OBJECT_INTERSECTIONS);

					for(size_t intersection = 0 ; intersection < (size_t)n_intersections; intersection++) {
						if(closest_intersection.x > intersections[intersection].x) {
							closest_intersection_object = obj;
							closest_intersection = intersections[intersection];
						}
					}
				}
			}
			if(closest_intersection_object < 0) continue;

			V3 intersection_normal;
			switch(scene[closest_intersection_object].type) {
				case SPHERE: 
					intersection_normal = normalize(vector_sub(closest_intersection, scene[closest_intersection_object].position)); 
					break;
				case PLANE:  
					intersection_normal = scene[closest_intersection_object].plane.normal; 
					break;
				default: 
					assert(0 && "Unreachable");
			}

			if(dot_product(intersection_normal, normalize(ray.direction)) > 0.0) {
				// the surface is facing away
				continue;
			}

			// find intersection -> light source vector
			Ray ray1 = {
				.direction = vector_sub(light_position, closest_intersection),
				.origin = closest_intersection,
			};

			// check for occlusions (shadows)
			int occluded = 0;
			for(size_t obj = 0; obj < scene_size; obj++) {
				if(obj != (size_t)closest_intersection_object && (scene[obj].type == SPHERE || scene[obj].type == PLANE)) {
					V3 intersections[MAX_OBJECT_INTERSECTIONS];
					int n_intersections = ray_object_intersect(scene[obj], ray1, intersections);
					assert(n_intersections >= 0);
					assert(n_intersections <= MAX_OBJECT_INTERSECTIONS);

					for(size_t intersection = 0 ; intersection < (size_t)n_intersections; intersection++) {
						if(point_in_segment(intersections[intersection], closest_intersection, light_position)) {
							// it is an occlusion
							occluded = 1;
							break;
						}
					}
				}
				if(occluded) break;
			}
			if(occluded) continue;
		
			float reflectance = dot_product(intersection_normal, normalize(ray1.direction));
			// assert(reflectance <= 1.0);
			
			if(reflectance < 0.0) {
				// the surface is facing away from the light source
				continue;
			}

			// assert(light.intensity <= 1.0 && light.intensity >= 0.0);
			Pixel p = PIXEL_RGB(
				(Pixel) (PIXEL_R(light.colour) * reflectance * light.intensity),
				(Pixel) (PIXEL_G(light.colour) * reflectance * light.intensity),
				(Pixel) (PIXEL_B(light.colour) * reflectance * light.intensity)
			);

			out_buffer[i * width + j] = p;
		}
	}
	

	return 0;
}

#define RENDERING_THREADS 16
#define SCENE_MAX_OBJECTS 1024

int objects_in_scene = 0;
Object scene[SCENE_MAX_OBJECTS];

sem_t sem_render_thread_start;
sem_t sem_render_thread_end;
void* new_render_thread(void* arg) {
	int i = *(int*) arg;
	assert(i >= 0);

	printf("thread %d created\n", i);

	while(1) {
		sem_wait(&sem_render_thread_start);
		sem_wait(&sem_render_thread_end);

		render_scene_threaded(scene, objects_in_scene, pixel_buffer, WIDTH, HEIGHT, RENDERING_THREADS, i);

		sem_post(&sem_render_thread_end);
	}

	return NULL;
}

int main() {
	// ** X display and window **
    Display *display;
    Window window;
	int screen;
	
    display = XOpenDisplay(NULL);
	if(!display) {
		printf("error: could not open X11 display\n");
		return 1;
	}

    screen = DefaultScreen(display);

    window = XCreateSimpleWindow(
		display, RootWindow(display, screen),
	   	0, 0, 
		WIDTH, HEIGHT, 
		WINDOW_BORDER_WIDTH,
		BlackPixel(display, screen), WhitePixel(display, screen)
	);

	// ** X events **
    XSelectInput(display, window, ExposureMask | KeyPressMask);
	

	// ** X graphics context **
	GC gc = XCreateGC(display, window, 0, NULL);
	

	// ** X image **
	XImage* image = XCreateImage(
		display, 
		DefaultVisual(display, screen),
	   	PIXEL_DEPTH, 
		ZPixmap, 
		0,
        (char*) pixel_buffer, 
		WIDTH, HEIGHT, 
		8 * sizeof(Pixel), 
		0
	);
	if(!image) {
		printf("error: could not create XImage.\n");
		return 1;
	}


	// create a scene
	scene[objects_in_scene++] = (Object){
		.type = CAMERA,
		.position = (V3){0, 0, 0},
		.camera = (Camera){
			.near_plane_distance = 2,
			.near_plane_width  = 1.8, // (float)WIDTH  / 1000.0,
			.near_plane_height = 1.8,  // (float)HEIGHT / 1000.0
			.far_plane_distance = 1000,
		},
	};
	scene[objects_in_scene++] = (Object){
		.type = LIGHT_SOURCE,
		.position  = {-10.0, 0.5, -20.0},
		.light_source = (LightSource){
			.colour = PIXEL_RGB(0x65, 0x32, 0x26),
			.intensity = 1.0
		},
	};
	scene[objects_in_scene++] = (Object){
		.type = PLANE,
		.position = {0.0, 0.0, 5.0},
		.plane = (Plane){
			.normal = {0.0, 0.0, -1.0},
		},
	};
	scene[objects_in_scene++] = (Object){
		.type = SPHERE,
		.position = {20.0, 0.0, 0.0},
		.sphere = (Sphere){
			.radius = 3.0,
		},
	};
	scene[objects_in_scene++] = (Object){
		.type = PLANE,
		.position = {30.0, 0.0, 0.0},
		.plane = (Plane){
			.normal = {-1.0, 0.0, 0.0},
		},
	};
	scene[objects_in_scene++] = (Object){
		.type = PLANE,
		.position = {0.0, 0.0, -25.0},
		.plane = (Plane){
			.normal = {0.0, 0.0, 1.0},
		},
	};
	scene[objects_in_scene++] = (Object){
		.type = PLANE,
		.position = {0.0, -7.0, 0.0},
		.plane = (Plane){
			.normal = {0.0, 1.0, 0.0},
		},
	};
	scene[objects_in_scene++] = (Object){
		.type = PLANE,
		.position = {0.0, 7.0, 0.0},
		.plane = (Plane){
			.normal = {0.0, -1.0, 0.0},
		},
	};

	// ** Rendering threads **
	sem_init(&sem_render_thread_start, 0, 0);
	sem_init(&sem_render_thread_end,   0, RENDERING_THREADS);

	pthread_t rendering_threads[RENDERING_THREADS];
	int rendering_threads_ids[RENDERING_THREADS];
	for(int i = 0; i < RENDERING_THREADS; i++) {
		rendering_threads_ids[i] = i;

		assert(!pthread_create(
			&rendering_threads[i],
			NULL,
			new_render_thread,
			(void*) &rendering_threads_ids[i]
		));
	}

	// ** Event loop **
    XMapWindow(display, window);

	// register interest in the delete window message
	Atom wm_delete_window = XInternAtom(display, "WM_DELETE_WINDOW", False);
	XSetWMProtocols(display, window, &wm_delete_window, 1);

	// find x11 events file descriptor
	struct pollfd xfd = {
		.fd = ConnectionNumber(display),
		.events = POLLIN,
	};

	// time setup
	struct timespec frame_end, frame_start;
	double dt, game_time;

	clock_gettime(CLOCK_MONOTONIC, &frame_start);

	int should_quit = 0;
	while (!should_quit) {
		{	
			while(XPending(display)) {
				XEvent event;
				XNextEvent(display, &event);

				switch (event.type) {
					case ClientMessage:
						if ((Atom)event.xclient.data.l[0] == wm_delete_window)
							should_quit = 1;
						break;
				}
			}

			int num_ready = poll(&xfd, 1, 0);
			if(num_ready < 0) {
				perror("poll()");
				break;
			}

			scene[3].position.z =  (sin(game_time * 2.0 * M_PI * 0.2)) * (scene[2].position.z - scene[3].sphere.radius);
			// render_scene(scene, objects_in_scene, pixel_buffer, WIDTH, HEIGHT);

			// start all the render threads
			for(int i = 0; i < RENDERING_THREADS; i++) {
				sem_post(&sem_render_thread_start);
			}

			// busy wait for all the render threads to end
			// TODO: figure out how to not busy wait
			int ended;
			do {
				sem_getvalue(&sem_render_thread_end, &ended);
				assert(ended >= 0);
			} while(ended < RENDERING_THREADS);


			XPutImage(display, window, gc, image, 0, 0, 0, 0, WIDTH, HEIGHT);
			XFlush(display);
		}

		clock_gettime(CLOCK_MONOTONIC, &frame_end);
		dt = (frame_end.tv_sec - frame_start.tv_sec) + (frame_end.tv_nsec - frame_start.tv_nsec) / 1e9;

		#ifdef MAX_FPS
		{
			double time_to_sleep = (double)1.0 / (double)MAX_FPS - dt;
			struct timespec ts;
			if(time_to_sleep > 0) {
				ts.tv_sec  = (time_t)time_to_sleep;
				ts.tv_nsec = (long)((time_to_sleep - ts.tv_sec) * 1e9);
				nanosleep(&ts, NULL);

				// get time after sleep
				clock_gettime(CLOCK_MONOTONIC, &frame_end);
				dt = (frame_end.tv_sec - frame_start.tv_sec) + (frame_end.tv_nsec - frame_start.tv_nsec) / 1e9;
			}
			
		}
		#endif
			
		frame_start = frame_end;
		game_time += dt;
		printf("FPS : %f\n", (double)1.0 / dt);
	}


    XDestroyWindow(display, window);
    XCloseDisplay(display);
    return 0;
}
