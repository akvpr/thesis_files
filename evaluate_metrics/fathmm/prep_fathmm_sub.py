if __name__ == '__main__':


	var_dict_pathogenic = {}
	var_dict_neutral = {}

	with open("../Pipeline/Data/p2_clean/testVariants_p2", 'r') as in_file:

		for line in in_file:

			line = line.strip()
			fields = line.split('\t')

			identifier = fields[0]
			sub = fields[1]
			label = fields[2]

			if label == 'neutral':
				var_dict = var_dict_neutral
			else:
				var_dict = var_dict_pathogenic

			if identifier not in var_dict:
				var_dict[identifier] = []

			var_dict[identifier].append(sub)


	with open('fathmm_submission_neutral', 'w') as out_file:

		for identifier in var_dict_neutral:

			out_string = identifier + '\t'
			var_string = ','.join(var_dict_neutral[identifier])
			out_string += var_string

			out_file.write(out_string + '\n')

	with open('fathmm_submission_pathogenic', 'w') as out_file:

		for identifier in var_dict_pathogenic:

			out_string = identifier + '\t'
			var_string = ','.join(var_dict_pathogenic[identifier])
			out_string += var_string

			out_file.write(out_string + '\n')