use crate::beam_utils::{PatientBox, TissueBox};
use rand::Rng;

#[derive(Debug, Clone, Copy)]
pub struct Vector {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

impl Vector {
    pub fn new(x_in: f32, y_in: f32, z_in: f32) -> Vector {
        Vector {
            x: x_in,
            y: y_in,
            z: z_in,
        }
    }

    pub fn dot(&self, v2: &Vector) -> f32 {
        self.x * v2.x + self.y * v2.y + self.z * v2.z
    }

    pub fn calculate_offset(&mut self, v2: &Vector) {
        self.x = self.x - v2.x as f32;
        self.y = self.y - v2.y as f32;
        self.z = self.z - v2.z as f32;
    }

    pub fn dist_to_beam(&self) -> f32 {
        let dist: f32 = (self.x.powf(2.0) + self.y.powf(2.0) + self.z.powf(2.0)) as f32;
        dist.sqrt()
    }

    pub fn dist_to_vector(&self, v2: &Vector) -> f32 {
        let dist: f32 = ((v2.x - self.x).powf(2.0)
            + (v2.y - self.y).powf(2.0)
            + (v2.z - self.z).powf(2.0)) as f32;
        dist.sqrt()
    }

    pub fn mult_vec(&mut self, val: f32) -> Vector {
        Vector {
            x: self.x * val,
            y: self.y * val,
            z: self.z * val,
        }
    }

    pub fn beam_direction(&self, tumour: &TissueBox) -> Vector {
        Vector {
            x: (tumour.x as f32 - self.x),
            y: (tumour.y as f32 - self.y),
            z: (tumour.z as f32 - self.z),
        }
    }

    pub fn crossover(&self, p2: &Vector, alpha: f32) -> Vector {
        Vector {
            x: crossover_val(&self.x, &p2.x, alpha),
            y: crossover_val(&self.y, &p2.y, alpha),
            z: crossover_val(&self.z, &p2.z, alpha),
        }
    }

    pub fn mutate(&mut self, mutation_bound: f32, patient: &PatientBox) {
        self.x = mutate_val(&self.x, mutation_bound, patient.x_size as f32);
        self.y = mutate_val(&self.y, mutation_bound, patient.y_size as f32);
        self.z = mutate_val(&self.z, mutation_bound, patient.z_size as f32);
    }
}

fn mutate_val(val: &f32, max_bound: f32, upper_bound: f32) -> f32 {
    if *val != 0.0 {
        let mut rng = rand::rng();
        let draw: f32 = rng.random_range(-max_bound..max_bound);
        let nv = (val + draw).max(0.0).min(upper_bound);
        nv
    } else {
        0.0
    }
}

fn crossover_val(x1: &f32, x2: &f32, alpha: f32) -> f32 {
    let mut nx: f32 = 0.0;
    if *x1 != 0.0 && *x2 != 0.0 {
        nx = alpha * x1 + (1.0 - alpha) * x2;
    }
    nx
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_crossover_val() {
        let ans = crossover_val(&3.0, &6.25, 0.5);
        assert_eq!(ans, 4.625);
        let ans = crossover_val(&0.0, &6.25, 0.5);
        assert_eq!(ans, 0.0);
    }
}
