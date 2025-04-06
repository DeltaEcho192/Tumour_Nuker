use crate::beam_utils::PatientBox;
use crate::beam_utils::generate_beam_entries;
use crate::vector::Vector;
use log::debug;
use rand::Rng;

#[derive(Clone)]
pub struct Indv {
    pub beams: Vec<Vector>,
    pub fitness: f32,
}
pub fn create_initial_population(pop_size: usize, patient_box: &PatientBox) -> Vec<Indv> {
    let mut population: Vec<Indv> = vec![];
    for _n in 0..pop_size {
        population.push(Indv {
            beams: generate_beam_entries(&patient_box),
            fitness: 0.0,
        });
    }
    population
}

pub fn selection(population: &Vec<Indv>, tournament_size: usize) -> Vec<Indv> {
    let mut rng = rand::rng();
    let mut selection: Vec<Indv> = vec![];
    let required_parents = population.len();
    if required_parents < tournament_size {
        panic!("Tournament size must be smaller than the amount of required parents");
    }
    let mut count = 0usize;
    loop {
        if count >= required_parents {
            break;
        }

        let mut max_indv: Option<&Indv> = None;
        let mut prev: Vec<usize> = vec![];
        let mut count_tour = 0usize;
        loop {
            if count_tour >= tournament_size {
                break;
            }
            let rng_idx = rng.random_range(0..required_parents);
            if prev.contains(&rng_idx) {
                continue;
            }

            if max_indv
                .unwrap_or(&Indv {
                    beams: vec![],
                    fitness: 0.0,
                })
                .fitness
                < population[rng_idx].fitness
            {
                debug!("New Fitness: {}", population[rng_idx].fitness);
                debug!(
                    "OLD Fitness: {}",
                    max_indv
                        .unwrap_or(&Indv {
                            beams: vec![],
                            fitness: 0.0,
                        })
                        .fitness
                );
                max_indv = Some(&population[rng_idx]);
            }
            count_tour += 1;
            prev.push(rng_idx);
        }

        selection.push(max_indv.unwrap().clone());
        count += 1;
    }

    selection
}

pub fn crossover(parent1: &Indv, parent2: &Indv) -> (Indv, Indv) {
    let mut rng = rand::rng();
    let alpha: f32 = rng.random_range(0.0..1.0);
    let child1 = Indv {
        beams: calculate_beam_crossover(&parent1.beams, &parent2.beams, alpha),
        fitness: 0.0,
    };
    let child2 = Indv {
        beams: calculate_beam_crossover(&parent2.beams, &parent1.beams, alpha),
        fitness: 0.0,
    };
    (child1, child2)
}

pub fn calculate_beam_crossover(
    p1_beams: &Vec<Vector>,
    p2_beams: &Vec<Vector>,
    alpha: f32,
) -> Vec<Vector> {
    if p1_beams.len() != p2_beams.len() {
        panic!("Amount of beams for crossover must be equal");
    }
    let mut new_beams: Vec<Vector> = vec![];
    for i in 0..p1_beams.len() {
        new_beams.push(p1_beams[i].crossover(&p2_beams[i], alpha));
    }
    new_beams
}
pub fn mutation() {
    todo!();
}

pub fn ga() {
    todo!();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_selection() {
        let population: Vec<Indv> = [
            Indv {
                beams: vec![],
                fitness: 20.0,
            },
            Indv {
                beams: vec![],
                fitness: 10.0,
            },
            Indv {
                beams: vec![],
                fitness: 50.0,
            },
        ]
        .to_vec();
        let ans = selection(&population, 2);
        assert_eq!(ans.len(), 3);
        for indv in ans {
            assert_ne!(indv.fitness, 10.0);
        }
    }
}
