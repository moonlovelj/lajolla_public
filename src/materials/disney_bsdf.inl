#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    (void)reflect; // silence unuse warning, remove this when implementing hw

    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specularTransmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    if (dot(vertex.geometric_normal, dir_in) <= 0)
    {
        return eval(DisneyGlass{ bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta }, dir_in, dir_out, vertex, texture_pool, dir)
            * (1.0 - metallic) * specularTransmission;
    }
    else
    {
        // Real specularTransmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
        // Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real specularTint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real sheenTint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real clearcoatGloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real eta = bsdf.eta;

        // Real diffuseWeight = (1 - metallic) * (1 - specularTransmission);
        // Real metalWeight = 1 - specularTransmission * (1 - metallic);
        // Real glassWeight = (1 - metallic) * specularTransmission;
        // Real clearcoatWeight = 0.25 * clearcoat;

        return eval(DisneyDiffuse{ bsdf.base_color, bsdf.roughness, bsdf.subsurface }, dir_in, dir_out, vertex, texture_pool, dir) * (1 - specularTransmission) * (1 - metallic) +
            eval(DisneySheen{ bsdf.base_color, bsdf.sheen_tint }, dir_in, dir_out, vertex, texture_pool, dir) * (1 - metallic) * sheen +
            eval(DisneyMetal{ bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.specular, bsdf.metallic, bsdf.specular_tint, bsdf.eta }, dir_in, dir_out, vertex, texture_pool, dir) * (1 - specularTransmission * (1 - metallic)) +
            eval(DisneyClearcoat{ bsdf.clearcoat_gloss }, dir_in, dir_out, vertex, texture_pool, dir) * 0.25 * clearcoat +
            eval(DisneyGlass{ bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta }, dir_in, dir_out, vertex, texture_pool, dir) * (1 - metallic) * specularTransmission;
    }
}

Real pdf_sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    (void)reflect; // silence unuse warning, remove this when implementing hw

    if (dot(vertex.geometric_normal, dir_in) <= 0)
    {
        return pdf_sample_bsdf(DisneyGlass{ bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta },
            dir_in, dir_out, vertex, texture_pool, dir);
    }
    else
    {
        Real specularTransmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real specularTint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real sheenTint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real clearcoatGloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real eta = bsdf.eta;

        Real diffuseWeight = (1 - metallic) * (1 - specularTransmission);
        Real metalWeight = 1 - specularTransmission * (1 - metallic);
        Real glassWeight = (1 - metallic) * specularTransmission;
        Real clearcoatWeight = 0.25 * clearcoat;
        Real numWeightInverse = 1 / (diffuseWeight + metalWeight + glassWeight + clearcoatWeight);
        diffuseWeight *= numWeightInverse;
        metalWeight *= numWeightInverse;
        glassWeight *= numWeightInverse;
        clearcoatWeight *= numWeightInverse;

        return pdf_sample_bsdf(DisneyDiffuse{ bsdf.base_color, bsdf.roughness, bsdf.subsurface }, dir_in, dir_out, vertex, texture_pool, dir) * diffuseWeight +
            pdf_sample_bsdf(DisneyMetal{ bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.specular, bsdf.metallic, bsdf.specular_tint, bsdf.eta }, dir_in, dir_out, vertex, texture_pool, dir) * metalWeight +
            pdf_sample_bsdf(DisneyClearcoat{ bsdf.clearcoat_gloss }, dir_in, dir_out, vertex, texture_pool, dir) * clearcoatWeight +
            pdf_sample_bsdf(DisneyGlass{ bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta }, dir_in, dir_out, vertex, texture_pool, dir) * glassWeight;
    }
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    if (dot(vertex.geometric_normal, dir_in) <= 0)
    {
        return sample_bsdf(DisneyGlass{ bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta }, 
            dir_in, vertex, texture_pool, rnd_param_uv, rnd_param_w, dir);
    }
    else
    {
        Real specularTransmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real specularTint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real sheenTint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real clearcoatGloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real eta = bsdf.eta;

        Real diffuseWeight = (1 - metallic) * (1 - specularTransmission);
        Real metalWeight = 1 - specularTransmission * (1 - metallic);
        Real glassWeight = (1 - metallic) * specularTransmission;
        Real clearcoatWeight = 0.25 * clearcoat;

        Real numWeightInverse = 1 / (diffuseWeight + metalWeight + glassWeight + clearcoatWeight);
        diffuseWeight *= numWeightInverse;
        metalWeight *= numWeightInverse;
        glassWeight *= numWeightInverse;
        clearcoatWeight *= numWeightInverse;

        if (rnd_param_w <= diffuseWeight)
        {
            return sample_bsdf(DisneyDiffuse{ bsdf.base_color, bsdf.roughness, bsdf.subsurface },
                dir_in, vertex, texture_pool, rnd_param_uv, rnd_param_w / diffuseWeight, dir);
        }
        else if (rnd_param_w <= diffuseWeight + metalWeight)
        {
            return sample_bsdf(DisneyMetal{ bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.specular, bsdf.metallic, bsdf.specular_tint, bsdf.eta },
                dir_in, vertex, texture_pool, rnd_param_uv, (rnd_param_w- diffuseWeight) / metalWeight, dir);
        }
        else if (rnd_param_w <= diffuseWeight + metalWeight + glassWeight)
        {
            return sample_bsdf(DisneyGlass{ bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta },
                dir_in, vertex, texture_pool, rnd_param_uv, (rnd_param_w - diffuseWeight - metalWeight) / glassWeight, dir);
        }
        else
        {
            return sample_bsdf(DisneyClearcoat{ bsdf.clearcoat_gloss},
                dir_in, vertex, texture_pool, rnd_param_uv, (rnd_param_w - diffuseWeight - metalWeight - glassWeight) / clearcoatWeight, dir);
        }
    }
}

TextureSpectrum get_texture_op::operator()(const DisneyBSDF &bsdf) const {
    return bsdf.base_color;
}
