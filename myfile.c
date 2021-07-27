#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
	float r[3], v[3];
} Particle;

const float METRE = 1.0;
const float SECOND = 1.0;
const float KG = 1.0;

const float B_MAN_CONST = 1.38e-23 * METRE * METRE * KG / (SECOND * SECOND);
const float MASS = 1.67e-27 * KG;
const float TEMPERATURE = 300;

const size_t PARTICLE_NUM = 1000;
const unsigned int SIM_PTOP_COLL_COUNT = 10000;
const float RADIUS = 0.01 * METRE;
const float DIAMETER = 2 * RADIUS;
const float PLANE_NEG_COORD = -5 * METRE;
const float PLANE_POS_COORD = 5 * METRE;
const float PLANE_POS[3] = {PLANE_POS_COORD, PLANE_POS_COORD, PLANE_POS_COORD};
const float PLANE_NEG[3] = {PLANE_NEG_COORD, PLANE_NEG_COORD, PLANE_NEG_COORD};

const float BIN_SIZE = 10 * METRE / SECOND;

float dot_product(float a[], float b[]) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

float sq_magnitude(float a[]) {
	return dot_product(a, a);
}

float discriminant(float a, float b, float c) {
	return b * b - 4 * a * c;
}

void add(float a[], float b[], float res[]) {
	res[0] = a[0] + b[0];
	res[1] = a[1] + b[1];
	res[2] = a[2] + b[2];
}

void subtract(float a[], float b[], float res[]) {
	res[0] = a[0] - b[0];
	res[1] = a[1] - b[1];
	res[2] = a[2] - b[2];
}

void resolve_parallel(float a[], float b[], float res[]) {
	float tmp = dot_product(a, b) / sq_magnitude(b);
	res[0] = tmp * b[0];
	res[1] = tmp * b[1];
	res[2] = tmp * b[2];
}

void initialize(Particle p[], const size_t n);
void simulate(Particle p[], const size_t n);
int print_data(Particle p[], const size_t n);

int main() {
	Particle p[PARTICLE_NUM];

	initialize(p, sizeof p / sizeof p[0]);
	simulate(p, sizeof p / sizeof p[0]);
	print_data(p, sizeof p / sizeof p[0]);

	return EXIT_SUCCESS;
}

void initialize(Particle p[], const size_t n) {
	size_t i, j;
	float v_rms = sqrt(3 * B_MAN_CONST * TEMPERATURE / MASS);
	float factor = 0;
	float temp[3];
	
	for(i = 0; i < n; i++) {
		for(j = 0; j < 3; j++) {
			p[i].v[j] = rand();
			factor += p[i].v[j] * p[i].v[j];
		}
		
		do {
			for (j = 0; j < 3; j++) {
				p[i].r[j] = (float) rand() / RAND_MAX * (PLANE_POS[j] - PLANE_NEG[j] - DIAMETER) + PLANE_NEG[j] + RADIUS;
			}
			for (j = 0; j < i; j++) {
				subtract(p[i].r, p[j].r, temp);
				if (sq_magnitude(temp) < DIAMETER * DIAMETER) break;
			}
		} while (j != i);
	}
	
	factor = v_rms * sqrt(n / factor);
	
	for(i = 0; i < n; i++) {
		for(j = 0; j < 3; j++) {
			p[i].v[j] *= factor;
		}
	}
}

void simulate(Particle p[], const size_t n) {
	size_t i, j, collpa, collpb;
	unsigned int ptop_coll_count = 0;
	float elapsed = 0, min_step, collision_time;
	int flag;
	float P, D;
	float r_ij[3], v_ij[3];
	float v_apar[3], v_bpar[3];

	while (ptop_coll_count < SIM_PTOP_COLL_COUNT) {
		/* Particles can either collide with the wall or with each other.
		 * Find out after how much time the next collision of any type
		 * takes place.
		 * Simultaneously figure out what particles undergo collision
		 * after said amount of time passes.
		 *
		 * TODO: Use Graph data structure to store all the colliding
		 * particles.
		*/
		
		min_step = INFINITY;
		for (i = 0; i < n; i++) {
			for (j = i + 1; j < n; j++) {
				subtract(p[i].r, p[j].r, r_ij);
				subtract(p[i].v, p[j].v, v_ij);
				P = dot_product(r_ij, v_ij);
				D = discriminant(
					sq_magnitude(v_ij),
					2 *(dot_product(r_ij, v_ij)),
					sq_magnitude(r_ij) - (DIAMETER * DIAMETER);
				if (P < 0 && D > 0) {
					collision_time = (-2 * P - sqrt(D))/(2 * sq_magnitude(v_ij));
					if (collision_time < min_step) {
						min_step = collision_time;
						collpa = i;
						collpb = j;
						flag = -1;
					}
				}
            }
			for (j = 0; j < 3; j++) {
				if (p[i].v[j] > 0) {
					collision_time = (PLANE_POS[j] - p[i].r[j] - RADIUS) / p[i].v[j];
					if (collision_time < min_step) {
						min_step = collision_time;
						collpa = i;
						flag = j;
					}
				} else if (p[i].v[j] < 0) {
					collision_time = (PLANE_NEG[j] - p[i].r[j] + RADIUS) / p[i].v[j];
					if (collision_time < min_step) {
						min_step = collision_time;
						collpa = i;
						flag = j;
					}
				}
			}
		}

		/* Now, move the particles to their corresponding locations they will
		 * move to, after min_step amount of time passes.
		*/
		for (i = 0; i < n; i++) {
			for (j = 0; j < 3; j++) {
				p[i].r[j] += p[i].v[j] * min_step;
			}
		}

		/* Update elapsed time */
		elapsed += min_step;

		/* Update the velocities of colliding particles to post-collision
		 * velocities.
		*/
		if (flag == 0) {
			p[collpa].v[0] = -p[collpa].v[0];
		} else if (flag == 1) {
			p[collpa].v[1] = -p[collpa].v[1];
		} else if (flag == 2) {
			p[collpa].v[2] = -p[collpa].v[2];
		} else {
			subtract(p[collpa].r, p[collpb].r, r_ij);
			
			resolve_parallel(p[collpa].v,  r_ij, v_apar);
			resolve_parallel(p[collpb].v, r_ij, v_bpar);
			
			/* v_ij is used as relative velocity to save memory */
			subtract(v_bpar, v_apar, v_ij);
			add(v_ij, p[collpa].v, p[collpa].v);
			subtract(p[collpb].v, v_ij, p[collpb].v);
			
			ptop_coll_count++;
			
			/* v1= (u2.r12)r12+ u1-(u1.r21)r21
			 * v2=u2-(u2.r12)r12+ (u1.r21)r21
			*/
			printf("Particle to Particle collision count: %u\t", ptop_coll_count);
			printf("Elapsed time: %f\n", elapsed);
		}
	}
	
}

int print_data(Particle p[], const size_t n) {
	for (int i = 0; i < n; i++) {
		printf("%f\n", sqrt(sq_magnitude(p[i].v)));
	}
}

