import argparse
import niio

parser = argparse.ArgumentParser()

parser.add_argument('--vector',help='Input label vector.',required=True,type=str)
parser.add_argument('--output',help='Output label Gifti file.',required=True,type=str)
parser.add_argument('--hemisphere',help='Left or Right hemisphere.',
	required=True,choices=['L','R'],type=str)

args = parser.parse_args()

hemimap = {'L': 'CortexLeft',
			'R': 'CortexRight'}

try:
	vector = niio.load(args.vector)
except:
	raise('File does not exist.')
else:
	niio.save(vector,args.output,hemimap[args.hemisphere])