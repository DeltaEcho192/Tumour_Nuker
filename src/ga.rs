use crate::beam_utils::PatientBox;
use crate::beam_utils::generate_beam_entries;
use crate::vector::Vector;
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
                println!("New Fitness: {}", population[rng_idx].fitness);
                println!("OLD Fitness: {}", max_indv
                .unwrap_or(&Indv {
                    beams: vec![],
                    fitness: 0.0,
                })
                .fitness);
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

pub fn crossover() {}

pub fn mutation() {}

pub fn ga() {}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_selection() {
        let population: Vec<Indv> = [Indv {beams: vec![], fitness: 20.0},Indv {beams: vec![], fitness: 10.0},Indv {beams: vec![], fitness: 50.0},].to_vec();  
        let ans = selection(&population, 2);
        assert_eq!(ans.len(), 3);
        for indv in ans {
            assert_ne!(indv.fitness, 10.0);
        }
    }
}
