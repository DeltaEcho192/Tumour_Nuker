use Performance_test::beam_utils;
use Performance_test::vector::Vector;
use std::thread;

fn loop_test() {
    println!("Run Loop test");
    const grid_x: i32 = 200;
    const grid_y: i32 = 100;
    const grid_z: i32 = 80;
    const total_size: usize = (grid_x * grid_y * grid_z) as usize;

    let mut inbeam_count: i64 = 0;
    let mut operationcounter: i64 = 0;
    let mut handle_vec = Vec::new();
    let beam_amt: i32 = 12;
    for i in 1..beam_amt {
        handle_vec.push(thread::spawn(|| {
            let mut dose: Vec<f32> = Vec::with_capacity(total_size);
            let beam_vec: Vector = Vector::new(0.8, 1.43, 6.3);

            let beam_ent_x: i32 = 2;
            let beam_ent_y: i32 = 2;
            let beam_ent_z: i32 = 0;

            let self_dot = beam_vec.dot(&beam_vec);

            let mut inbeam_count: i64 = 0;
            let mut operationcounter: i64 = 0;
            for x in 0..grid_x {
                for y in 0..grid_y {
                    for z in 0..grid_z {
                        let mut dose_val: f32 = 0.0;
                        let mut vector = Vector::new(x as f32, y as f32, z as f32);
                        vector.calculate_offset(beam_ent_x, beam_ent_y, beam_ent_z);
                        let dist = vector.dist_to_beam();
                        let dot_prod = vector.dot(&beam_vec);
                        let projection_point = vector.mult_vec(dot_prod / self_dot);
                        let project_dist = vector.dist_to_vector(&projection_point);
                        if project_dist <= 0.75 {
                            inbeam_count += 1;
                            dose_val = project_dist;
                        }
                        //let index: usize = (x + grid_y * (y + grid_z * z)) as usize;
                        dose.push(dose_val);
                        operationcounter += 1;
                    }
                }
            }
        }));
    }

    for handle in handle_vec {
        handle.join().unwrap();
    }
}

fn main() {
    loop_test();
    println!("Hello, world!");
}
