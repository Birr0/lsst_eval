import requests
import pyvo
import os
from dotenv import load_dotenv
from tqdm import tqdm


load_dotenv()

def get_auth():
    ssotap_url = os.getenv("SSO_TAP_URL")
    TOKEN = os.getenv("TOKEN")

    s = requests.Session()
    s.headers["Authorization"] = "Bearer " + TOKEN
    auth = pyvo.auth.authsession.AuthSession()
    auth.credentials.set("lsst-token", s)

    auth.add_security_method_for_url(ssotap_url, "lsst-token")
    auth.add_security_method_for_url(ssotap_url + "/sync", "lsst-token")
    auth.add_security_method_for_url(ssotap_url + "/async", "lsst-token")
    auth.add_security_method_for_url(ssotap_url + "/tables", "lsst-token")
    return auth

def service():
    return pyvo.dal.TAPService(os.getenv("SSO_TAP_URL"), get_auth())

def query(stmt):
    return service().search(
        stmt
    )
def submit_job_list(statements):
    jobs = set()
    for stmt in statements:
        job = service().submit_job(
            stmt
        )
        job.run()
        jobs.add(job)
    return jobs

def check_jobs(jobs):
    completed = 0
    queue = 0

    for result in tqdm(list(jobs)):
        result.executionduration = 600
        if result.phase == "COMPLETED":
            completed += 1
        else:
            queue += 1
    print(f"Completed {completed} jobs. {queue} jobs still the queue.")
    return queue, jobs


def get_job_results(jobs):
    resp = []
    for result in tqdm(list(jobs)):
            result.executionduration = 600
            if result.phase == "COMPLETED":
                resp.append(result.fetch_result())
            else:
                # DEAL with other errors
                pass
    return resp

def cleanup_jobs():
    # To-Do
    return