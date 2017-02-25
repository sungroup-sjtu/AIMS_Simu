from flask import request
from ..models import *
from .. import db
from typing import Dict, List
from sqlalchemy import and_, or_


class ComputeAction():
    def init_from_json(self, compute_json) -> int:
        if compute_json['type'] == 'unary':
            type = Compute.Type.UNARY
        elif compute_json['type'] == 'binary':
            type = Compute.Type.BINARY
        else:
            raise Exception('Invalid compute type')

        try:
            simulations = self.get_simulation_types(compute_json['properties'])
        except Exception as e:
            raise Exception(str(e))

        if len(simulations) == 0 or len(compute_json['molecules']) == 0:
            raise Exception('Invalid information')

        compute = Compute()
        compute.web_id = compute_json['id']
        compute.web_user_id = compute_json['user_id']
        compute.web_ip = request.remote_addr
        compute.type = type

        db.session.add(compute)
        db.session.flush()

        if type == Compute.Type.UNARY:
            for molecule_json_dict in compute_json['molecules']:
                for simulation in simulations:
                    task = ComputeUnary()
                    task.compute_id = compute.id
                    task.web_molecule_id = molecule_json_dict['id']
                    task.web_molecule_smiles = molecule_json_dict['smiles']
                    task.simulation_id = simulation.id
                    task.t_min = compute_json['t_min']
                    task.t_max = compute_json['t_max']
                    task.p_min = compute_json['p_min']
                    task.p_max = compute_json['p_max']
                    db.session.add(task)

        elif type == Compute.Type.BINARY:
            for molecule_json_dict in compute_json['molecules']:
                for molecule2_json_dict in compute_json['molecule2s']:
                    for simulation in simulations:
                        task = ComputeBinary()
                        task.compute_id = compute.id
                        task.web_molecule_id = molecule_json_dict['id']
                        task.web_molecule_smiles = molecule_json_dict['smiles']
                        task.web_molecule2_id = molecule2_json_dict['id']
                        task.web_molecule2_smiles = molecule2_json_dict['smiles']
                        task.simulation_id = simulation.id
                        task.t_min = compute_json['t_min']
                        task.t_max = compute_json['t_max']
                        task.p_min = compute_json['p_min']
                        task.p_max = compute_json['p_max']
                        db.session.add(task)

        try:
            db.session.commit()
        except Exception as e:
            db.session.rollback()
            raise Exception(str(e))
        else:
            return compute.id

    def get_simulation_types(self, properties_json_dict: List[Dict]) -> List[Simulation]:
        simulations = []
        for property_dict in properties_json_dict:
            name = property_dict['name']
            abbr = property_dict['abbr']

            property = Property.query.join(Simulation).filter(
                or_(
                    Property.name.ilike(name),
                    Property.abbr.ilike(abbr),
                )
            ).first()

            if property == None:
                raise Exception('Invalid property: ' + name)

            if property.simulation not in simulations:
                simulations.append(property.simulation)
        return simulations
