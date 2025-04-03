use crate::beam_utils::{PatientBox, TissueBox};
use std::cmp;

pub struct MaskHolder {
    pub masks: Vec<Mask>,
}

pub struct Mask {
    pub x0: i64,
    pub x1: i64,
    pub y0: i64,
    pub y1: i64,
    pub z0: i64,
    pub z1: i64,
}

impl Mask {
    pub fn from_tissue_box(t_box: &TissueBox, p_box: &PatientBox) -> Mask {
        Mask {
            x0: cmp::max(0, t_box.x - t_box.x_width / 2),
            x1: cmp::min(p_box.x_size, t_box.x + t_box.x_width / 2),
            y0: cmp::max(0, t_box.y - t_box.y_width / 2),
            y1: cmp::min(p_box.y_size, t_box.y + t_box.y_width / 2),
            z0: cmp::max(0, t_box.z - t_box.z_width / 2),
            z1: cmp::min(p_box.z_size, t_box.z + t_box.z_width / 2),
        }
    }

    pub fn bound_check(&self, x: i64, y: i64, z: i64) -> bool {
        if x >= self.x0
            && x <= self.x1
            && y >= self.y0
            && y <= self.y1
            && z >= self.z0
            && z <= self.z1
        {
            true
        } else {
            false
        }
    }
}
