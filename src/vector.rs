#[derive(Debug, Clone)]
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
}
