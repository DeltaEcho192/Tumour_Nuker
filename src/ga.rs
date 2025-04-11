use crate::beam_utils::PatientBox;
use crate::beam_utils::TissueBox;
use crate::beam_utils::{ComputeDoseParamsIter, compute_costIter, compute_dose_iter, generate_beam_entries};
use crate::mask::Mask;
use crate::vector::Vector;
use std::sync::{RwLock, Arc};
use log::debug;
use rand::Rng;

#[derive(Clone)]
pub struct Indv {
    pub beams: Vec<Vector>,
    pub fitness: f32,
}

impl Indv {
    pub fn calculate_fitness<const N_SIZE: usize>(
        &mut self,
        patient: &PatientBox,
        tumour: &TissueBox,
        mask_holder: &Vec<Mask>,
    ) {
        let mut dose_params: ComputeDoseParamsIter<{ N_SIZE }> = ComputeDoseParamsIter {
            patient_box: patient.clone(),
            beams: self.beams.clone(),
            tumour: tumour.clone(),
            dose_matrix: vec![0f32; N_SIZE].try_into().unwrap(),
        };
        compute_dose_iter(&mut dose_params);
        self.fitness = compute_costIter(&mut dose_params, &mask_holder);
    }
    pub fn mutation(&mut self, patient: &PatientBox, mutation_prop: f32, mutation_bound: f32) {
        let mut rng = rand::rng();
        let draw: f32 = rng.random_range(0.0..1.0);
        if draw <= mutation_prop {
            for beam in &mut self.beams {
                beam.mutate(mutation_bound, &patient);
            }
        }
    }
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

        let mut min_indv: Option<&Indv> = None;
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

            if min_indv
                .unwrap_or(&Indv {
                    beams: vec![],
                    fitness: 100000000.0,
                })
                .fitness
                > population[rng_idx].fitness
            {
                debug!("New Fitness: {}", population[rng_idx].fitness);
                debug!(
                    "OLD Fitness: {}",
                    min_indv
                        .unwrap_or(&Indv {
                            beams: vec![],
                            fitness: 0.0,
                        })
                        .fitness
                );
                min_indv = Some(&population[rng_idx]);
            }
            count_tour += 1;
            prev.push(rng_idx);
        }

        selection.push(min_indv.unwrap().clone());
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

const MUTATION_PROB: f32 = 0.25;
const MUTATION_BOUND: f32 = 10.0;

pub fn ga<const N_SIZE: usize>(
    population_size: usize,
    generations: usize,
    patient: PatientBox,
    tumour: TissueBox,
    mask_holder: Vec<Mask>,
    tournament_size: usize,
) {
    let mut population = create_initial_population(population_size, &patient);
    let mut best_in_gen: Vec<Indv> = vec![];
    for generation in 0..generations {
        for indv in &mut population {
            indv.calculate_fitness::<{ N_SIZE }>(&patient, &tumour, &mask_holder);
        }
        let reproduce_pop = selection(&population, tournament_size);
        let mut new_pop: Vec<Indv> = vec![];
        let max_idx: usize = reproduce_pop.len();
        let mut idx: usize = 0;
        loop {
            if (idx + 1) >= max_idx {
                break;
            }
            let (mut child1, mut child2) = crossover(&reproduce_pop[idx], &reproduce_pop[idx + 1]);
            child1.mutation(&patient, MUTATION_PROB, MUTATION_BOUND);
            child2.mutation(&patient, MUTATION_PROB, MUTATION_BOUND);
            new_pop.push(child1);
            new_pop.push(child2);
            idx += 2;
        }
        let mut min_indv: Option<&Indv> = None;
        for indv in &population {
            if min_indv
                .unwrap_or(&Indv {
                    beams: vec![],
                    fitness: 10000000000.0,
                })
                .fitness
                > indv.fitness
            {
                min_indv = Some(&indv);
            }
        }
        new_pop[0] = min_indv.unwrap().clone();
        best_in_gen.push(min_indv.unwrap().clone());
        println!("Best Solution in Generation {} : {}", generation, min_indv.unwrap().fitness);
        population = new_pop;
        println!("Population Size: {}", population.len());
        println!(
            "Percentage Done: {} %",
            (generation as f32 / generations as f32) * 100.0
        );
    }

    let mut min_indv: Option<&Indv> = None;
    for indv in &best_in_gen {
        if min_indv
            .unwrap_or(&Indv {
                beams: vec![],
                fitness: 100000000.0,
            })
        .fitness
        > indv.fitness
        {
            min_indv = Some(&indv);
        }
    }

    println!("Best Solution Len: {}", best_in_gen.len());
    println!("Best Solution: {}", min_indv.unwrap().fitness);
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
            assert_ne!(indv.fitness, 50.0);
        }
    }
}
