import requests
import asyncio

async def get_gene_name(gene_id: str)->dict:
	'''
	Queries wormbase for a given gene id.
	Return json containing gene name.
	'''
	loop = asyncio.get_event_loop()
	headers = {"Accept": "application/json"}

	uri = f"http://www.wormbase.org//rest/field/gene/{gene_id}/name"
	res = await loop.run_in_executor(None,
		partial(requests.get,
			uri,
			headers=headers))
	if res:
		return res.json()

if __name__ == "__main__":
	gene_id = "WBGene00002915"
	response = asyncio.run(get_gene_name(gene_id))
	print(response)