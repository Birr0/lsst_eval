import requests
import pyvo
import os
from dotenv import load_dotenv

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