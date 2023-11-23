#pragma once

// The simplest volumetric renderer: 
// single absorption only homogeneous volume
// only handle directly visible light sources
Spectrum vol_path_tracing_1(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);

    RayDifferential ray_diff;
    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
    if (!vertex_) {
        return make_zero_spectrum();
    }

    PathVertex vertex = *vertex_;

    Spectrum Le = make_zero_spectrum();
    if (is_light(scene.shapes[vertex.shape_id])) {
        Le = emission(vertex, -ray.dir, scene);
    }

    Spectrum transmittance = make_const_spectrum(1);
    if (vertex.exterior_medium_id >= 0)
    {
        transmittance = exp(-get_sigma_a(scene.media[vertex.exterior_medium_id], vertex.position) * distance(vertex.position, ray.org));
    }

    return Le * transmittance;
}

// The second simplest volumetric renderer: 
// single monochromatic homogeneous volume with single scattering,
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_2(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);

    RayDifferential ray_diff;
    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
    
    Real t_hit = infinity<Real>();
    if (vertex_) {
        t_hit = distance((*vertex_).position, ray.org);
    }
    
    if (scene.camera.medium_id < 0)
    {
        Spectrum Le = make_zero_spectrum();
        if (vertex_ && is_light(scene.shapes[(*vertex_).shape_id])) {
            Le = emission(*vertex_, -ray.dir, scene);
        }
        return Le;
    }
    
    Spectrum sigma_t = get_sigma_a(scene.media[scene.camera.medium_id], Vector3{0, 0, 0}) + get_sigma_s(scene.media[scene.camera.medium_id], Vector3{0, 0, 0});
    Real u = next_pcg32_real<Real>(rng);
    Real t = -log(1-u) / sigma_t.x;
    if (t < t_hit)
    {
        Spectrum trans_pdf = exp(-sigma_t * t) * sigma_t;
        Spectrum transmittance = exp(-sigma_t * t);
        Vector3 p = ray.org + t * ray.dir;
        Spectrum sigma_s = get_sigma_s(scene.media[scene.camera.medium_id], p);
        Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
        Real light_w = next_pcg32_real<Real>(rng);
        Real shape_w = next_pcg32_real<Real>(rng);
        int light_id = sample_light(scene, light_w);
        const Light &light = scene.lights[light_id];
        PointAndNormal point_on_light =
            sample_point_on_light(light, p, light_uv, shape_w, scene);

        Real G = 0;
        Vector3 dir_light = normalize(point_on_light.position - p);
        Ray shadow_ray{p, dir_light, 
                               get_shadow_epsilon(scene),
                               (1 - get_shadow_epsilon(scene)) *
                                   distance(point_on_light.position, p)};
        if (!occluded(scene, shadow_ray)) {
            // geometry term is cosine at v_{i+1} divided by distance squared
            // this can be derived by the infinitesimal area of a surface projected on
            // a unit sphere -- it's the Jacobian between the area measure and the solid angle
            // measure.
            G = max(-dot(dir_light, point_on_light.normal), Real(0)) /
                distance_squared(point_on_light.position, p);
        }

        Real p1 = light_pmf(scene, light_id) *
                pdf_point_on_light(light, point_on_light, p, scene);

        if (G > 0 && p1 > 0)
        {
            Spectrum L = emission(light, -dir_light, Real(0), point_on_light, scene);
            PhaseFunction phase_function = get_phase_function(scene.media[scene.camera.medium_id]);
            Spectrum phase_eval = eval(phase_function,  -ray.dir, dir_light);
            return (transmittance / trans_pdf) * sigma_s * exp(-sigma_t * distance(point_on_light.position, p)) * L / p1 * G * phase_eval;
        }
        
        return make_zero_spectrum();
    }
    else
    {
        Spectrum Le = make_zero_spectrum();
        if (is_light(scene.shapes[(*vertex_).shape_id])) {
            Le = emission((*vertex_), -ray.dir, scene);
        }
        Spectrum trans_pdf = exp(-sigma_t * t_hit);
        Spectrum transmittance = exp(-sigma_t * t_hit);
        return transmittance / trans_pdf * Le;
    }
}

// The third volumetric renderer (not so simple anymore): 
// multiple monochromatic homogeneous volumes with multiple scattering
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_3(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    int current_medium_id = scene.camera.medium_id;
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;
    Spectrum current_path_throughput = make_const_spectrum(1);

    while (true)
    {
        bool scatter = false;
        std::optional<PathVertex> vertex_ = intersect(scene, ray);
        Spectrum transmittance = make_const_spectrum(1);
        Spectrum trans_pdf = make_const_spectrum(1);
        Real t_hit = infinity<Real>();
        if (vertex_) {
            t_hit = distance((*vertex_).position, ray.org);
        }
        if (current_medium_id > -1)
        {
            Spectrum sigma_t = get_sigma_a(scene.media[current_medium_id], Vector3{0, 0, 0}) + get_sigma_s(scene.media[current_medium_id], Vector3{0, 0, 0});
            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1-u) / sigma_t.x;
            if (t < t_hit)
            {
                scatter = true;
                trans_pdf = exp(-sigma_t * t) * sigma_t;
                transmittance = exp(-sigma_t * t);
            }
            else
            {
                trans_pdf = exp(-sigma_t * t_hit);
                transmittance = exp(-sigma_t * t_hit);
            }

            current_path_throughput *= (transmittance / trans_pdf);
            ray.org = ray.org + t * ray.dir;
        }

        if (!scatter)
        {
            Spectrum Le = make_zero_spectrum();
            if (vertex_ && is_light(scene.shapes[(*vertex_).shape_id])) {
                Le = emission((*vertex_), -ray.dir, scene);
            }
            radiance += current_path_throughput * Le;
        }

        if (bounces == scene.options.max_depth - 1 && scene.options.max_depth != -1)
        {
            break;
        }

        if (!scatter && vertex_ && (*vertex_).material_id == -1)
        {
            if ( (*vertex_).interior_medium_id !=  (*vertex_).exterior_medium_id)
            {
                if (dot(ray.dir, (*vertex_).geometric_normal) > 0)
                {
                    current_medium_id = (*vertex_).exterior_medium_id;
                }
                else
                {
                    current_medium_id = (*vertex_).interior_medium_id;
                }
            }
            bounces += 1;
            continue;
        }

        if (scatter)
        {
            PhaseFunction phase_function = get_phase_function(scene.media[current_medium_id]);
            std::optional<Vector3> next_dir = sample_phase_function(phase_function, -ray.dir, Vector2 (next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)));
            current_path_throughput *=eval(phase_function, -ray.dir, *next_dir) / pdf_sample_phase(phase_function, -ray.dir, *next_dir) * get_sigma_s(scene.media[current_medium_id], Vector3{0, 0, 0});
            ray.dir = *next_dir;
        }
        else
        {
            break;
        }

        Real rr_prob = 1;
        if (bounces >= scene.options.rr_depth)
        {
            rr_prob = min(current_path_throughput.x, 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob)
            {
                break;
            }
            else
            {
                current_path_throughput /= rr_prob;
            }
        }

        bounces += 1;
    }
    
    return radiance;
}

// The fourth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// still no surface lighting
Spectrum vol_path_tracing_4(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}

// The fifth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing_5(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}

// The final volumetric renderer: 
// multiple chromatic heterogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing(const Scene &scene,
                          int x, int y, /* pixel coordinates */
                          pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}
