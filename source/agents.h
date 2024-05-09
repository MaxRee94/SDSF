#pragma once
#include "helpers.h"
#include "grid_agent.forward.h"
#include "kernel.h"


using namespace help;

class Strategy {
public:
	Strategy() = default;
	Strategy(string _vector, float _seed_mass, float _diaspore_mass, int _no_seeds_per_diaspore, float _seed_tspeed,
		float _pulp_to_seed_ratio, float _recruitment_probability
	) :
		vector(_vector), seed_mass(_seed_mass), diaspore_mass(_diaspore_mass), no_seeds_per_diaspore(_no_seeds_per_diaspore),
		seed_tspeed(_seed_tspeed), pulp_to_seed_ratio(_pulp_to_seed_ratio), recruitment_probability(_recruitment_probability)
	{}
	void print() {
		printf("id: %d, seed_mass: %f, diaspore_mass: %f, no_seeds_per_diaspore: %d, vector: %s, pulp to seed ratio: %f, seed terminal speed: %f, germination prob: %f\n",
			id, seed_mass, diaspore_mass, no_seeds_per_diaspore, vector.c_str(), pulp_to_seed_ratio, seed_tspeed, recruitment_probability
		);
	}
	int id = -1;
	float seed_mass = 0;
	float diaspore_mass = 0;
	int no_seeds_per_diaspore = 0;
	float seed_tspeed = 0;
	float pulp_to_seed_ratio = 0;
	float recruitment_probability = 0;
	string vector = "none";
};


class StrategyGenerator {
public:
	StrategyGenerator() = default;
	StrategyGenerator(map<string, map<string, float>> user_parameters) {
		for (auto& [trait, value_range] : user_parameters) {
			ProbModel prob_model;
			string distribution_type = distribution_types[value_range["distribution_type"]];
			if (distribution_type == "uniform") {
				prob_model = ProbModel(value_range["min"], value_range["max"]);
			}
			else if (distribution_type == "linear") {
				prob_model = ProbModel(value_range["q1"], value_range["q2"], value_range["min"], value_range["max"]);
			}
			else if (distribution_type == "normal") {
				prob_model = ProbModel(value_range["mean"], value_range["stdev"], 0);
			}
			else if (distribution_type == "discrete") {
				prob_model = ProbModel(value_range["probability0"], value_range["probability1"], value_range["probability2"]);
			}
			trait_distributions[trait] = prob_model;
		}
	}
	string pick_vector() {
		int vector = trait_distributions["vector"].sample();
		if (vector == 0) return "linear";
		else if (vector == 1) return "wind";
		else return "animal";
	}
	float compute_wing_mass(float cumulative_seed_mass) {
		// Estimate wing mass based on correlation we found in seed- and wing measurement data taken from (Greene and Johnson, 1993)
		return max(0, (cumulative_seed_mass - 0.03387f) / 3.6609f);
	}
	float compute_diaspore_mass(float no_seeds_per_diaspore, float seed_mass, string vector, float fruit_pulp_mass, float &pulp_to_seed_ratio) {
		float cumulative_seed_mass = seed_mass * no_seeds_per_diaspore;
		float diaspore_mass;
		if (vector == "wind") {
			float wing_mass = compute_wing_mass(cumulative_seed_mass);
			diaspore_mass = cumulative_seed_mass + wing_mass;
			//printf("cumulative seed mass: %f, wing mass: %f \n", cumulative_seed_mass, wing_mass);
		}
		else if (vector == "animal") {
			diaspore_mass = cumulative_seed_mass + fruit_pulp_mass;
		}
		else {
			diaspore_mass = cumulative_seed_mass;
		}
		pulp_to_seed_ratio = fruit_pulp_mass / cumulative_seed_mass;
		return diaspore_mass;
	}
	float calculate_tspeed(float diaspore_mass) {
		// Compute terminal descent velocity based on correlation presented by (Greene and Johnson, 1993)
		return 0.501f * pow(diaspore_mass * 1000, 0.174);
	}
	int sample_no_seeds_per_diaspore() {
		int no_seeds_per_diaspore = round(trait_distributions["no_seeds_per_diaspore"].sample());
		if (no_seeds_per_diaspore < 1) no_seeds_per_diaspore = 1; // Ensure at least one seed per diaspore
		return no_seeds_per_diaspore;
	}
	float sample_seed_mass(string vector) {
		float seed_mass = 1.0f;
		if (vector == "wind") seed_mass = trait_distributions["seed_mass_wind"].sample();
		if (vector == "animal") seed_mass = trait_distributions["seed_mass_animal"].sample();
		return seed_mass;
	}
	float sample_fruit_pulp_mass() {
		return trait_distributions["fruit_pulp_mass"].sample();
	}
	float calculate_recruitment_probability(float seed_mass) {
		return 0.0385 * log(seed_mass) + 0.224; // Fitted to data from Barczyk et al (2024), see file 'seed weight vs seedling success.xlsx'
	}
	void generate(Strategy &strategy) {
		int no_seeds_per_diaspore = sample_no_seeds_per_diaspore();
		string vector = pick_vector();
		float seed_mass = sample_seed_mass(vector);
		float fruit_pulp_mass = sample_fruit_pulp_mass();
		float pulp_to_seed_ratio;
		float diaspore_mass = compute_diaspore_mass(no_seeds_per_diaspore, seed_mass, vector, fruit_pulp_mass, pulp_to_seed_ratio);
		float seed_tspeed = calculate_tspeed(diaspore_mass);
		float recruitment_probability = calculate_recruitment_probability(seed_mass);
		strategy = Strategy(vector, seed_mass, diaspore_mass, no_seeds_per_diaspore, seed_tspeed, pulp_to_seed_ratio, recruitment_probability);
	}
	void mutate(Strategy& strategy, float mutation_rate) {
		bool do_mutation = help::get_rand_float(0, 1) < mutation_rate;
		if (do_mutation) {
			int trait_idx = help::get_rand_int(0, 3);
			float fruit_pulp_mass = 0;
			if (trait_idx == 0) {
				strategy.seed_mass = sample_seed_mass(strategy.vector);
			}
			else if (trait_idx == 1) {
				strategy.no_seeds_per_diaspore = sample_no_seeds_per_diaspore();
			}
			else if (trait_idx == 2) {
				fruit_pulp_mass = sample_fruit_pulp_mass();
			}
			else if (trait_idx == 3) {
				strategy.vector = pick_vector();
			}
			float pulp_to_seed_ratio;
			strategy.diaspore_mass = compute_diaspore_mass(strategy.no_seeds_per_diaspore, strategy.seed_mass, strategy.vector, fruit_pulp_mass, pulp_to_seed_ratio);
			strategy.seed_tspeed = calculate_tspeed(strategy.diaspore_mass);
			strategy.pulp_to_seed_ratio = pulp_to_seed_ratio;
		}
	}
	map<string, ProbModel> trait_distributions;
	map<int, string> distribution_types = { { 0, "uniform" }, {1, "linear"}, {2, "normal"}, {3, "discrete"} };
};


class Tree {
public:
	Tree() = default;
	Tree(int _id, pair<float, float> _position, float _radius, float _radius_tmin1, float& seed_bearing_threshold) : 
		position(_position), radius(_radius), radius_tmin1(_radius_tmin1) 
	{
		id = _id;
		init();
	};
	Tree(int _id, pair<float, float> _position, float _radius, float _radius_tmin1, int _life_phase, float &seed_bearing_threshold) :
		position(_position), radius(_radius), radius_tmin1(_radius_tmin1), life_phase(_life_phase)
	{
		id = _id;
		init();
		//printf("Tree created with id: %i, radius: %f, radius_tmin1: %f, stem dbh: %f, bark thickness: %f, LAI: %f \n", id, radius, radius_tmin1, stem_dbh, bark_thickness, LAI);
	};
	bool operator==(const Tree& tree) const
	{
		return id == tree.id;
	}
	void init() {
		stem_dbh = get_stem_dbh_from_radius();
		bark_thickness = get_bark_thickness();
		LAI = get_LAI();
		//printf("LAI: %f, stem dbh: %f, radius: %f, radius_tmin1: %f \n", LAI, stem_dbh, radius, radius_tmin1);
	}
	bool is_within_radius(pair<float, float> pos2) {
		float dist = help::get_dist(position, pos2);
		return dist < radius;
	}
	int get_life_phase(float& seed_bearing_threshold) {
		if (radius > seed_bearing_threshold) return 2;
		else return 0;
	}
	float get_bark_thickness() {
		return 0.31 * pow(stem_dbh, 1.276); // From Hoffmann et al (2012), figure 5a. Bark thickness in mm.
	}
	float get_stem_dbh_from_radius() {
		return pow(10.0f, (log10(radius + radius) + 0.12) / 0.63); // From Antin et al (2013), figure 1, topright panel (reordered equation). Stem dbh in cm.
	}
	float get_stem_dbh_increment() {
		// TODO: Ensure increment is influenced by seed size.
		float stem_increment = 3.0f * (1.0f - exp(-0.118 * stem_dbh));	// I = I_max(1-e^(-g * D)), from Hoffman et al (2012), Appendix 1, page 1.
																		// Current stem dbh and increment in cm.

		stem_increment *= (5.0f - LAI) * 0.2f;	// LAI-dependent growth reduction to introduce density dependence. Multiplication factor: ((LAI_max - LAI) / LAI_max) 
												// From Hoffman et al (2012), Appendix 2.

		return stem_increment;
	}
	float get_survival_probability(float& fire_resistance_argmin, float& fire_resistance_argmax, float& fire_resistance_stretch) {
		return help::get_sigmoid(bark_thickness, fire_resistance_argmin, fire_resistance_argmax, fire_resistance_stretch); // Sigmoid based on Hoffman et al (2012), figure 2a.
	}
	float get_leaf_area() {
		return 0.147 * pow(stem_dbh, 2.053); // From Hoffman et al (2012), figure 5b. Leaf area in m^2
	}
	float get_ground_area() {
		return M_PI * (radius * radius);
	}
	float get_LAI() {
		float ground_area = get_ground_area();
		float leaf_area = get_leaf_area();
		return leaf_area / ground_area;
	}
	bool survives_fire(float &fire_resistance_argmin, float &fire_resistance_argmax, float &fire_resistance_stretch) {
		float survival_probability = get_survival_probability(fire_resistance_argmin, fire_resistance_argmax, fire_resistance_stretch);
		return help::get_rand_float(0.0f, 1.0f) > survival_probability;
	}
	float get_crown_radius_from_dbh() {
		float diameter = pow(10.0f, -0.12 + 0.63 * log10(stem_dbh)); // From Antin et al (2013), figure 1, topright panel. Stem dbh in cm.
		return diameter * 0.5f;
	}
	float get_height() {
		float ln_dbh = log(stem_dbh);
		return exp(0.865 + 0.760 * ln_dbh - 0.0340 * (ln_dbh * ln_dbh)); // From Chave et al (2014), equation 6a. Value of E obtained by calculating mean of E
																		 // values for bistable study sites.
	}
	void grow(float &seed_bearing_threshold) {
		age++;
		radius_tmin1 = radius;
		stem_dbh += get_stem_dbh_increment();
		radius = get_crown_radius_from_dbh();
		life_phase = get_life_phase(seed_bearing_threshold);
		bark_thickness = get_bark_thickness();
		LAI = get_LAI();
		height = get_height();
	}
	float radius = -1;
	float radius_tmin1 = -1;
	float stem_dbh = 0.1;
	float bark_thickness = 0;
	float LAI = 0;
	float height = 0;
	pair<float, float> position = pair(0, 0);
	int id = -1;
	int age = -1;
	int life_phase = 0;
	int last_mortality_check = 0;
};


class Crop {
public:
	Crop() = default;
	Crop(Strategy &_strategy, Tree& tree, float _mass_budget_factor) {
		strategy = _strategy;
		seed_mass = _strategy.seed_mass;
		mass_budget_factor = _mass_budget_factor;
		origin = tree.position;
		id = tree.id;
	}
	void compute_no_seeds(Tree &tree) {
		no_seeds = 49324 * exp(0.0491 * tree.stem_dbh);	// Fitted to averages generated by model: #seedlings = STR' * (dbh / 30)^alpha, modified from Uriarte et al(2005),
														// where STR' = STR / (influence of seed mass on STR (fitted to normalized STR and seed mass values)), so that
														// the model is seed mass-normalized).
		no_seeds *= 0.0771 * exp(-4.956 * (strategy.seed_mass / 3.56)); // Re-introduce (normalized) seed mass influence.
		no_seeds *= strategy.recruitment_probability;
	}
	void compute_no_diaspora() {
		no_diaspora = no_seeds / strategy.no_seeds_per_diaspore;
		no_seeds = no_diaspora * strategy.no_seeds_per_diaspore;
	}
	void update(Tree& tree) {
		if (tree.id != id) printf("\n\n--------------- CROP ID (%i) DOES NOT MATCH TREE ID (%i) ---------------\n\n", id, tree.id);
		compute_no_seeds(tree);
		compute_no_diaspora();
	}
	int no_seeds = 0;
	int no_diaspora = 0;
	float seed_mass = 0;
	float mass_budget_factor = 0;
	Strategy strategy;
	pair<float, float> origin = pair<float, float>(0, 0);
	int id = -1;
};


class Population {
public:
	Population() = default;
	Population(float _max_radius, float _cellsize, float _radius_q1, float _radius_q2, float _mass_budget_factor, 
		map<string, map<string, float>> strategy_parameters, float _mutation_rate, float _seed_bearing_threshold
	) : max_radius(_max_radius), cellsize(_cellsize), radius_q1(_radius_q1), radius_q2(_radius_q2),
		mass_budget_factor(_mass_budget_factor), seed_bearing_threshold(_seed_bearing_threshold)
	{
		strategy_generator = StrategyGenerator(strategy_parameters);
		radius_probability_model = help::LinearProbabilityModel(radius_q1, radius_q2, 0, max_radius);
		mutation_rate = _mutation_rate;
	}
	Tree* add(pair<float, float> position, Strategy* _strategy = 0, float radius = -2) {
		// Create tree
		if (radius == -1) radius = max_radius;
		else if (radius == -2) radius = radius_probability_model.linear_sample();
		float radius_tmin1 = radius * 0.9; // TEMP. TODO: Make dependent on growth curve.
		Tree tree(no_created_trees + 1, position, radius, radius_tmin1, 1, seed_bearing_threshold);
		members[tree.id] = tree;
		no_created_trees++;

		// Create strategy
		Strategy strategy;
		if (_strategy != nullptr) {
			strategy = *_strategy;
			strategy_generator.mutate(strategy, mutation_rate);
		}
		else {
			strategy_generator.generate(strategy);
		} 
		strategy.id = tree.id;

		// Create crop
		Crop crop(strategy, tree, mass_budget_factor);
		crops[tree.id] = crop;

		// Create custom kernel
		Kernel* _kernel = get_kernel(strategy.vector);
		Kernel kernel;
		if (strategy.vector == "wind") {
			kernel = Kernel(
				tree.id, _kernel->dist_max, _kernel->wspeed_gmean, _kernel->wspeed_stdev, _kernel->wind_direction,
				_kernel->wind_direction_stdev, strategy.seed_tspeed, max_radius * 4
			);
		}
		kernels_individual[tree.id] = kernel;

		return &members[tree.id];
	}
	Kernel* add_kernel(string tree_dispersal_vector, Kernel &kernel) {
		kernels[tree_dispersal_vector] = kernel;
		return &kernel;
	}
	Tree* get(int id) {
		return &members[id];
	}
	void get(vector<int>& ids, vector<Tree*> &trees) {
		for (int id : ids) trees.push_back(get(id));
	}
	Crop* get_crop(int id) {
		return &crops[id];
	}
	Kernel* get_kernel(string vector) {
		return &kernels[vector];
	}
	Kernel* get_kernel(int id) {
		return &kernels_individual[id];
	}
	int size() {
		return members.size();
	}
	void remove(Tree* tree) {
		remove(tree->id);
	}
	void remove(int id) {
		members.erase(id);
		crops.erase(id);
		kernels_individual.erase(id);
	}
	bool is_population_member(Tree* tree) {
		auto it = members.find(tree->id);
		return (it != members.end());
	}
	bool is_population_member(int tree_id) {
		return members.find(tree_id) != members.end();
	}
	unordered_map<int, Tree> members;
	unordered_map<int, Crop> crops;
	unordered_map<string, Kernel> kernels;
	unordered_map<int, Kernel> kernels_individual;
	help::LinearProbabilityModel radius_probability_model;
	StrategyGenerator strategy_generator;
	float max_radius = 0;
	float cellsize = 0;
	float radius_q1 = 0;
	float radius_q2 = 0;
	float seed_mass = 0;
	float mutation_rate = 0;
	Tree removed_tree;
	int no_created_trees = 0;
	float mass_budget_factor = 0;
	float seed_bearing_threshold = 0;
};
